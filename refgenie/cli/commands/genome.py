"""The `genome` command group: models and handlers."""

from collections.abc import Callable
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.framework import CliList
from refgenie.cli.commands.helpers import data_channel_hint, genome_digests_from
from refgenie.cli.errors import EXIT_INVALID_INPUT, EXIT_NOT_FOUND, fail
from refgenie.logger import logger
from refgenie.managers.collection_sync import collection_field, sync_source
from refgenie.managers.sources import make_source
from refgenie.managers.sources.genomes import RefgetStoreSource


class GenomeInitModel(BaseModel):
    """genome init: initialize a genome and build its fasta asset.

    Creates a genome from a FASTA file, server, or store, and builds its
    fasta asset (fa, fai, chrom.sizes) unless --no-build is specified. For a
    server or store genome, the build first downloads its sequences.
    """

    name: CliList = Field(
        description="One or more alias names for the genome.",
        validation_alias=AliasChoices("n", "name"),
    )
    fasta: Path | None = Field(
        default=None,
        description="Path to local FASTA file.",
    )
    server: str | None = Field(
        default=None,
        description="URL of a refgenie server to init from.",
    )
    store: str | None = Field(
        default=None,
        description="URL of a RefgetStore to init directly from (requires --digest or --namespace).",
    )
    namespace: str | None = Field(
        default=None,
        description="Namespace for alias lookup when using --store.",
    )
    digest: str | None = Field(
        default=None,
        description="Seqcol digest of the genome.",
    )
    description: str = Field(
        default="",
        description="Genome description (e.g. 'Human genome build 38').",
        validation_alias=AliasChoices("d", "description"),
    )
    species: str | None = Field(
        None,
        description="Species name (e.g. 'Homo sapiens').",
        validation_alias=AliasChoices("s", "species"),
    )
    fhr: Path | None = Field(
        default=None,
        description="Path to an FHR .fhr.json metadata file to apply after init. "
        "When given it is authoritative for description/species and writes the "
        "RefgetStore sidecar.",
    )
    force: bool = Field(
        False,
        description="Allow re-initialization of an existing genome (adds new aliases).",
        validation_alias=AliasChoices("f", "force"),
    )
    build: bool = Field(
        True,
        description="Automatically build the fasta asset after genome initialization. "
        "Only runs when a 'fasta' recipe is registered (via a synced data channel); "
        "otherwise a message is printed and the build is skipped. Use --no-build to "
        "skip the auto-build entirely.",
    )


class GenomeSetMetadataModel(BaseModel):
    """genome set-metadata: apply FHR metadata to an already-registered genome.

    Reads a normalized FHR .fhr.json file (a path token -- no shell-quoting of
    free-text metadata) and updates the genome's description/species columns and
    RefgetStore sidecar with no rebuild.
    """

    name: str | None = Field(
        default=None,
        description="Genome alias to update.",
        validation_alias=AliasChoices("n", "name"),
    )
    digest: str | None = Field(
        default=None,
        description="Genome seqcol digest to update (alternative to --name).",
    )
    fhr: Path = Field(
        description="Path to the FHR .fhr.json metadata file to apply.",
    )


class GenomeListModel(BaseModel):
    """genome list: list all genomes."""

    pass


class GenomeRemoveModel(BaseModel):
    """genome remove: remove a genome and all its assets."""

    genome: CliList | None = Field(
        default=None,
        description="Genome alias(es) to remove.",
    )
    genome_digest: CliList | None = Field(
        default=None,
        description="Genome digest(s) to remove, in place of aliases.",
        validation_alias=AliasChoices("genome-digest"),
    )
    force: bool = Field(
        False,
        description="Do not prompt before action.",
        validation_alias=AliasChoices("f", "force"),
    )


class GenomeBrowseModel(BaseModel):
    """genome browse: browse available genomes on a refgenie server or RefgetStore."""

    server_url: str | None = Field(
        default=None,
        description="URL of a refgenie server or a RefgetStore to browse. "
        "Defaults to subscribed server(s).",
    )
    page: int = Field(
        default=0,
        description="Page number for paginated results.",
    )
    page_size: int = Field(
        default=20,
        description="Number of results per page.",
        validation_alias=AliasChoices("page-size"),
    )


class GenomeSyncModel(BaseModel):
    """genome sync: bulk-register all genomes from a remote source."""

    server_url: str | None = Field(
        default=None,
        description="URL of a remote source (refgenie server or RefgetStore) to "
        "sync from. Defaults to all subscribed server(s).",
        validation_alias=AliasChoices("server-url"),
    )
    page_size: int = Field(
        default=1000,
        description="Number of collections to request per page.",
        validation_alias=AliasChoices("page-size"),
    )


class GenomeParser(BaseModel):
    """Intermediate parser for genome subcommands."""

    init: CliSubCommand[GenomeInitModel] = Field(
        description="Initialize a genome from a FASTA file or remote server."
    )
    set_metadata: CliSubCommand[GenomeSetMetadataModel] = Field(
        alias="set-metadata",
        description="Apply FHR metadata to an already-registered genome (no rebuild).",
    )
    list: CliSubCommand[GenomeListModel] = Field(description="List all genomes.")
    remove: CliSubCommand[GenomeRemoveModel] = Field(
        description="Remove a genome and all its assets."
    )
    browse: CliSubCommand[GenomeBrowseModel] = Field(
        description="Browse available genomes on a refgenie server or RefgetStore."
    )
    sync: CliSubCommand[GenomeSyncModel] = Field(
        description="Bulk-register all genomes from subscribed servers or a remote source."
    )


def _load_fhr_file(path: Path) -> dict:
    """Read an FHR .fhr.json file into a dict, failing cleanly if it is absent."""
    import json

    if not Path(path).exists():
        fail(f"FHR file not found: {path}")
    with open(path) as fh:
        return json.load(fh)


def _apply_fhr_file(refgenie, digest: str, path: Path) -> None:
    """Apply the FHR record at ``path`` to ``digest`` via GenomeManager.apply_fhr."""
    refgenie.genome.apply_fhr(digest, _load_fhr_file(path))
    logger.info(f"Applied FHR metadata from {path} to {digest}")


def handle_genome_init(cmd, refgenie) -> None:
    import tempfile

    from refgenie.exceptions import BuildFailedError

    if sum(s is not None for s in (cmd.fasta, cmd.server, cmd.store)) != 1:
        fail("Provide exactly one of --fasta, --server, or --store")

    with tempfile.TemporaryDirectory() as cache_dir:
        source, digest = None, cmd.digest
        if cmd.server:
            # make_source detects the backend, so a store URL passed to --server works too.
            try:
                source = make_source(cmd.server)
            except ConnectionError as e:
                fail(str(e))
        elif cmd.store:
            if not digest and not cmd.namespace:
                fail("--store requires either --digest or --namespace")
            source = RefgetStoreSource(cmd.store, Path(cache_dir))
            if digest is None:
                digest = source.resolve_alias_in_namespace(cmd.name[0], cmd.namespace)
                if digest is None:
                    fail(
                        f"Alias '{cmd.name[0]}' not found in namespace '{cmd.namespace}'",
                        EXIT_NOT_FOUND,
                    )
                logger.info(
                    f"Resolved alias '{cmd.name[0]}' in namespace '{cmd.namespace}' to digest: {digest}"
                )
        try:
            digest, created = refgenie.build.initialize_and_build(
                genome_names=cmd.name,
                fasta_file_path=cmd.fasta,
                source=source,
                digest=digest,
                description=cmd.description,
                species_name=cmd.species,
                use_existing=cmd.force,
                build_fasta=cmd.build,
            )
        except BuildFailedError as e:
            fail(f"{e} Retry with: refgenie build {cmd.name[0]}/fasta")

    # FHR overrides --description/--species, and is applied to existing genomes too.
    if cmd.fhr is not None:
        _apply_fhr_file(refgenie, digest, cmd.fhr)

    if not created:
        logger.info(f"Genome already exists: {digest}")
        return
    logger.info(f"Genome initialized: {digest} (aliases: {', '.join(cmd.name)})")
    if cmd.build and not refgenie.recipe.exists("fasta"):
        logger.info(data_channel_hint(f"refgenie build {cmd.name[0]}/fasta"))


def handle_genome_set_metadata(cmd, refgenie) -> None:
    from refgenie.exceptions import MissingAliasError

    if bool(cmd.digest) == bool(cmd.name):
        fail("Provide exactly one of --name (an alias) or --digest", EXIT_INVALID_INPUT)
    try:
        [digest] = genome_digests_from(
            refgenie, [cmd.name] if cmd.name else None, [cmd.digest] if cmd.digest else None
        )
    except MissingAliasError:
        fail(f"Genome '{cmd.name}' not found", EXIT_NOT_FOUND)
    _apply_fhr_file(refgenie, digest, cmd.fhr)


def handle_genome_list(cmd, refgenie) -> None:
    table = refgenie.genome.table()
    rprint(table)


def handle_genome_remove(cmd, refgenie) -> None:
    from rich.prompt import Confirm

    digests = genome_digests_from(refgenie, cmd.genome, cmd.genome_digest)
    if not digests:
        fail("Name the genome(s) to remove with --genome or --genome-digest", EXIT_INVALID_INPUT)
    for digest in digests:
        if not cmd.force:
            if not Confirm.ask(
                f"Are you sure you want to remove genome '{digest}' and all its assets?"
            ):
                logger.info("Aborted by user. Genome not removed.")
                continue
        refgenie.genome.remove(digest)
        logger.info(f"Removed genome: {digest}")


def handle_genome_browse(cmd, refgenie) -> None:
    if cmd.server_url:
        server_urls = [cmd.server_url]
    else:
        server_urls = refgenie.servers.subscriptions()
        if not server_urls:
            fail(
                "No server URL provided and no server subscriptions found.\n"
                "Either provide a URL: refgenie genome browse --server-url <URL>\n"
                "Or subscribe to a server: refgenie subscribe -s <URL>"
            )

    failures = 0
    for server_url in server_urls:
        try:
            source = make_source(server_url)
            result = source.list_collections(page=cmd.page, page_size=cmd.page_size)
        except ConnectionError as e:
            logger.error(f"Could not connect to server: {server_url}\n{e}")
            failures += 1
            continue

        if result is None:
            logger.error(f"No response from server: {server_url}")
            failures += 1
            continue

        items = result.get("results", [])
        if not items:
            if len(server_urls) > 1:
                print(f"No collections found on {server_url}.")
            else:
                print("No collections found.")
            continue

        if len(server_urls) > 1:
            print(f"\n--- {server_url} ---")

        for item in items:
            digest = collection_field(item, "digest") or "?"
            n_seqs = collection_field(item, "n_sequences")
            names_digest = collection_field(item, "names_digest") or ""
            print(f"{digest}\t{'?' if n_seqs is None else n_seqs} sequences\t{names_digest}")

    if failures:
        fail(f"{failures} of {len(server_urls)} server(s) could not be browsed")


def handle_genome_sync(cmd, refgenie) -> None:
    """Bulk-register all collections from subscribed servers or a remote source."""
    if cmd.server_url:
        server_urls = [cmd.server_url]
    else:
        server_urls = refgenie.servers.subscriptions()
        if not server_urls:
            fail(
                "No server URL provided and no server subscriptions found.\n"
                "Either provide a URL: refgenie genome sync --server-url <URL>\n"
                "Or subscribe to a server: refgenie subscribe -s <URL>"
            )

    total_registered = 0
    failures = 0
    for server_url in server_urls:
        logger.info(f"Syncing genomes from: {server_url}")
        try:
            source = make_source(server_url)
        except Exception as e:  # noqa: BLE001 - one unreachable server never stops a sync
            logger.error(f"Could not connect to {server_url}: {e}")
            failures += 1
            continue
        report = sync_source(refgenie, source, server_url, cmd.page_size)
        total_registered += report.registered
        failures += report.failures
        logger.info(f"  {report.registered} new genome(s) from {server_url}")

    logger.info(f"Synced {total_registered} new genome(s) total")
    if failures:
        fail(f"{failures} genome(s)/source(s) failed to sync")


GENOME_DISPATCH: dict[type, Callable] = {
    GenomeInitModel: handle_genome_init,
    GenomeSetMetadataModel: handle_genome_set_metadata,
    GenomeListModel: handle_genome_list,
    GenomeRemoveModel: handle_genome_remove,
    GenomeBrowseModel: handle_genome_browse,
    GenomeSyncModel: handle_genome_sync,
}


def handle_genome_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = GENOME_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown genome subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
