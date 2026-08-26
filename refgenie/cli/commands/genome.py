"""The `genome` command group: models and handlers."""

from collections.abc import Callable
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.commands.framework import CliList
from refgenie.cli.commands.helpers import data_channel_hint
from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.logger import logger
from refgenie.managers.sources import make_source
from refgenie.managers.sources.genomes import RefgetStoreSource


class GenomeInitModel(BaseModel):
    """genome init: initialize a genome and build its fasta asset.

    Creates a genome from a FASTA file, server, or store. When initialized
    from a FASTA file, automatically builds the fasta asset (fa, fai,
    chrom.sizes) unless --no-build is specified.
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

    genome: CliList = Field(
        description="Genome digest(s) or alias(es) to remove.",
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


def _make_source(url: str):
    """Create a RemoteGenomeSource from a URL, auto-detecting backend."""
    return make_source(url)


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
    sources = [s for s in (cmd.fasta, cmd.server, cmd.store) if s is not None]
    if len(sources) != 1:
        fail("Provide exactly one of --fasta, --server, or --store")

    if cmd.fasta:
        digest, created = refgenie.genome.initialize_genome(
            fasta_file_path=cmd.fasta,
            description=cmd.description,
            alias_names=cmd.name,
            species_name=cmd.species,
            use_existing=cmd.force,
        )
    elif cmd.server:
        # Route through make_source so --server and --store differ only in which
        # detection branch wins; a store URL passed to --server still works.
        try:
            source = _make_source(cmd.server)
        except ConnectionError as e:
            fail(str(e))
        digest, created = refgenie.genome.initialize_genome(
            source=source,
            digest=cmd.digest,
            description=cmd.description,
            alias_names=cmd.name,
            species_name=cmd.species,
            use_existing=cmd.force,
        )
    elif cmd.store:
        if not cmd.digest and not cmd.namespace:
            fail("--store requires either --digest or --namespace")
        import tempfile
        from pathlib import Path as _Path

        cache_dir = _Path(tempfile.mkdtemp())
        source = RefgetStoreSource(cmd.store, cache_dir)
        digest = cmd.digest
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
        digest, created = refgenie.genome.initialize_genome(
            source=source,
            digest=digest,
            description=cmd.description,
            alias_names=cmd.name,
            species_name=cmd.species,
            use_existing=cmd.force,
        )
    else:
        fail("Provide exactly one of --fasta, --server, or --store")

    # Apply FHR metadata if provided. Authoritative over --description/--species,
    # and applied on the use_existing/--force branch too -- the only way an
    # already-registered genome's metadata is refreshed without a rebuild.
    if getattr(cmd, "fhr", None) is not None:
        _apply_fhr_file(refgenie, digest, cmd.fhr)

    if created:
        logger.info(f"Genome initialized: {digest} (aliases: {', '.join(cmd.name)})")
        if cmd.build:
            if not refgenie.recipe.exists("fasta"):
                logger.info(
                    "Genome initialized. Skipping automatic fasta build: no 'fasta' recipe "
                    "is registered. "
                    + data_channel_hint(f"refgenie build {cmd.name[0]}/fasta")
                )
            else:
                # build_asset signals failure by RETURNING None, not by raising, so
                # the return value must be checked -- an except clause alone silently
                # reports success for a build that failed.
                try:
                    logger.info(f"Building fasta asset for {cmd.name[0]}...")
                    built = refgenie.build_asset(
                        recipe_name="fasta",
                        genome_name=cmd.name[0],
                        asset_group_name="fasta",
                    )
                except Exception as e:
                    # Recipes shell out to arbitrary tools, so a build can raise
                    # anything. Swallow it here only to reach the non-zero exit
                    # below with a legible message.
                    built = None
                    logger.error(f"Fasta asset build raised for {cmd.name[0]}: {e}")

                if built is None:
                    # Exit non-zero. A genome whose fasta asset did not build is not
                    # usable, and callers gate on this command's exit status.
                    fail(
                        f"Genome '{cmd.name[0]}' ({digest}) was registered but its fasta "
                        f"asset failed to build. Retry with: "
                        f"refgenie build {cmd.name[0]}/fasta"
                    )
                logger.info(f"Fasta asset built successfully for {cmd.name[0]}")
    else:
        logger.info(f"Genome already exists: {digest}")


def handle_genome_set_metadata(cmd, refgenie) -> None:
    from refgenie.exceptions import MissingAliasError

    if not cmd.digest and not cmd.name:
        fail("Provide --name or --digest")
    digest = cmd.digest
    if digest is None:
        try:
            digest = refgenie.alias.resolve(cmd.name)
        except MissingAliasError:
            fail(f"Genome alias '{cmd.name}' not found", EXIT_NOT_FOUND)
    _apply_fhr_file(refgenie, digest, cmd.fhr)


def handle_genome_list(cmd, refgenie) -> None:
    table = refgenie.genome.table()
    rprint(table)


def handle_genome_remove(cmd, refgenie) -> None:
    from rich.prompt import Confirm

    from refgenie.exceptions import MissingAliasError

    for genome_identifier in cmd.genome:
        try:
            digest = refgenie.alias.resolve(genome_identifier)
        except MissingAliasError:
            # Not an alias -- assume the user passed a digest directly.
            digest = genome_identifier
        if not cmd.force:
            if not Confirm.ask(
                f"Are you sure you want to remove genome '{digest}' and all its assets?"
            ):
                logger.info("Aborted by user. Genome not removed.")
                continue
        refgenie.genome.remove(digest)
        logger.info(f"Removed genome: {digest}")


def handle_genome_browse(cmd, refgenie) -> None:
    # Resolve server URLs: explicit argument or fall back to subscriptions
    if cmd.server_url:
        server_urls = [cmd.server_url]
    else:
        server_urls = refgenie.sources.get_subscriptions()
        if not server_urls:
            fail(
                "No server URL provided and no server subscriptions found.\n"
                "Either provide a URL: refgenie genome browse --server-url <URL>\n"
                "Or subscribe to a server: refgenie subscribe -s <URL>"
            )

    failures = 0
    for server_url in server_urls:
        try:
            source = _make_source(server_url)
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
            digest = (
                item.get("digest", "?") if isinstance(item, dict) else getattr(item, "digest", "?")
            )
            n_seqs = (
                item.get("n_sequences", "?")
                if isinstance(item, dict)
                else getattr(item, "n_sequences", "?")
            )
            names_digest = (
                item.get("names_digest", "")
                if isinstance(item, dict)
                else getattr(item, "names_digest", "")
            )
            print(f"{digest}\t{n_seqs} sequences\t{names_digest}")

    if failures:
        fail(f"{failures} of {len(server_urls)} server(s) could not be browsed")


# The one store alias namespace sync carries into refgenie's flat alias space.
# It is the only namespace that is 1:1 with collections (a curated, unique name
# per collection), so it maps cleanly onto flat aliases. The others (accession,
# insdc, refseq, genome_assembly, ...) hold coarse or duplicated names whose
# flat-space assignment would be order-dependent -- and genome_assembly names
# like 'hg38' would collide with genomes registered by builds.
SYNC_ALIAS_NAMESPACE = "name"


def _syncable_alias_names(refgenie, digest: str, alias_pairs: list[tuple[str, str]]) -> list[str]:
    """Filter a collection's store aliases down to names sync may register.

    Keeps only the SYNC_ALIAS_NAMESPACE namespace, and NEVER repoints: an alias
    that already resolves to a different genome (e.g. 'hg38' owned by a built
    genome) is skipped with a warning -- AliasManager.add would silently steal
    it otherwise. Aliases already pointing at this digest pass through, so
    repeat syncs can repair ones that went missing.
    """
    names = []
    for namespace, name in alias_pairs:
        if namespace != SYNC_ALIAS_NAMESPACE:
            continue
        try:
            existing = refgenie.alias.resolve(name)
        except Exception:  # noqa: BLE001 - unresolvable reads as unclaimed
            existing = None
        if existing is not None and existing != digest:
            logger.warning(
                f"Sync skipping alias '{name}' for {digest}: it already points "
                f"to {existing} and sync never repoints an alias."
            )
            continue
        if name not in names:
            names.append(name)
    return names


def handle_genome_sync(cmd, refgenie) -> None:
    """Bulk-register all collections from subscribed servers or a remote source.

    For each source, paginates through its collections and registers each as a
    genome carrying the store's curated aliases and FHR metadata. Registration
    is idempotent (use_existing=True): existing genomes are kept, with missing
    aliases repaired and FHR re-applied.
    """

    # Resolve source URLs: explicit argument or fall back to subscriptions.
    if cmd.server_url:
        server_urls = [cmd.server_url]
    else:
        server_urls = refgenie.sources.get_subscriptions()
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
            source = _make_source(server_url)
        except Exception as e:
            logger.error(f"Could not connect to {server_url}: {e}")
            failures += 1
            continue

        registered_here = 0
        page = 0
        while True:
            try:
                result = source.list_collections(page=page, page_size=cmd.page_size)
            except Exception as e:
                logger.error(f"Failed to list collections from {server_url}: {e}")
                failures += 1
                break
            if not result:
                break
            items = result.get("results", [])
            if not items:
                break

            for item in items:
                coll_digest = (
                    item.get("digest") if isinstance(item, dict) else getattr(item, "digest", None)
                )
                if not coll_digest:
                    continue
                description = (
                    item.get("description", "")
                    if isinstance(item, dict)
                    else getattr(item, "description", "")
                )
                try:
                    alias_pairs = source.get_collection_aliases(coll_digest)
                    alias_names = _syncable_alias_names(refgenie, coll_digest, alias_pairs)
                    digest, created = refgenie.genome.initialize_genome(
                        source=source,
                        digest=coll_digest,
                        description=description or "",
                        alias_names=alias_names,
                        use_existing=True,
                    )
                    if created:
                        registered_here += 1
                        logger.info(f"Registered genome: {digest}")
                except Exception as e:
                    logger.error(f"Failed to register {coll_digest}: {e}")
                    failures += 1
                    continue

                # apply_fhr is the only writer of the faceted metadata columns
                # (common_name/taxon_id/assembly_*). Applied on every sync, not
                # just creation, so store-side metadata edits propagate; a
                # collection with no FHR still registers with those columns
                # null. A failing apply counts as a sync failure but must not
                # be mistaken for a failed registration -- the genome row stands.
                try:
                    fhr = source.get_collection_fhr(coll_digest)
                    if fhr:
                        refgenie.genome.apply_fhr(digest, fhr)
                except Exception as e:
                    logger.error(f"Genome {digest} registered but FHR apply failed: {e}")
                    failures += 1

            if len(items) < cmd.page_size:
                break
            page += 1

        total_registered += registered_here
        logger.info(f"  {registered_here} new genome(s) from {server_url}")

    logger.info(f"Synced {total_registered} new genome(s) total")
    if failures:
        # A partial sync is a failure: callers gate downstream work on this
        # command's status, and exiting 0 hides genomes that never registered.
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
