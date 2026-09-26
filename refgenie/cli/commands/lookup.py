"""The commands that resolve one thing and print it -- `seek`, `seekr`, `id`,
`compare`: models and handlers.
"""

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg
from rich import print as rprint

from refgenie.cli.framework import CliList
from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _ASSET_REGISTRY_PATHS_SEEK_DESCRIPTION,
    _GENOME_DIGEST_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    _lookup_one,
    genome_digest_from,
    genome_ref,
    resolve_transient_servers,
)
from refgenie.cli.errors import EXIT_INVALID_INPUT, EXIT_NOT_FOUND, fail, run_over_paths
from refgenie.logger import logger


class SeekModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description=_ASSET_REGISTRY_PATHS_SEEK_DESCRIPTION,
    )
    check_exists: bool = Field(
        default=False,
        description="Whether the returned asset path should be checked for existence on disk.",
        validation_alias=AliasChoices("e", "check-exists"),
    )
    abs: bool = Field(
        default=False,
        description="Return the digest-addressed content path under data/ instead "
        "of the human-readable alias path.",
        validation_alias=AliasChoices("abs"),
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION,
        validation_alias=AliasChoices("genome-digest"),
    )


class SeekrModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description=_ASSET_REGISTRY_PATHS_SEEK_DESCRIPTION,
    )
    genome_server: CliList | None = Field(
        default=None,
        description=_GENOME_SERVER_TRANSIENT_DESCRIPTION,
        validation_alias=AliasChoices("s", "genome-server"),
    )
    append_server: bool = Field(
        default=False,
        description=_APPEND_SERVER_DESCRIPTION,
        validation_alias=AliasChoices("p", "append-server"),
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION,
        validation_alias=AliasChoices("genome-digest"),
    )


class IdModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        default_factory=list,
        description="One or more registry paths: genome alias (e.g. hg38) for "
        "genome digest, or asset path (e.g. hg38/fasta or hg38/fasta:default) "
        "for asset digest.",
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION
        + " Alone, prints the digest if the genome is known.",
        validation_alias=AliasChoices("genome-digest"),
    )
    verbose: bool = Field(
        default=False,
        description="Show detailed genome metadata (sequence count, total length, source).",
        validation_alias=AliasChoices("v", "verbose"),
    )
    validate_store: bool = Field(
        default=False,
        description="Verify the genome's RefgetStore exists and is valid.",
        validation_alias=AliasChoices("validate-store"),
    )
    remote: bool = Field(
        default=False,
        description="Query subscribed seqcolapi servers if the --genome-digest genome "
        "is not found locally.",
        validation_alias=AliasChoices("remote"),
    )
    info: bool = Field(
        default=False,
        description="Given a digest, show its aliases and metadata.",
        validation_alias=AliasChoices("info"),
    )


class CompareModel(BaseModel):
    genomes: CliPositionalArg[list[str]] = Field(
        default_factory=list,
        description="Genome aliases to compare. Together with --genome-digest, "
        "exactly two genomes.",
    )
    genome_digest: CliList | None = Field(
        default=None,
        description="Genome digests to compare, in place of aliases.",
        validation_alias=AliasChoices("genome-digest"),
    )


def handle_seek(cmd, refgenie) -> None:
    force_exists = getattr(cmd, "check_exists", False)
    abs_path = getattr(cmd, "abs", False)

    def lookup(parsed):
        if cmd.genome_digest:
            return refgenie.asset.seek(
                genome_ref(None, cmd.genome_digest),
                parsed.asset_group,
                asset_name=parsed.asset,
                seek_key_name=parsed.seek_key,
                force_exists=force_exists,
                abs_path=abs_path,
            )
        return refgenie.asset.seek_components(parsed, force_exists=force_exists, abs_path=abs_path)

    run_over_paths(
        cmd.asset_registry_paths,
        lambda raw: _lookup_one(refgenie, raw, lookup, cmd.genome_digest),
    )


def handle_seekr(cmd, refgenie) -> None:
    server_urls = resolve_transient_servers(cmd, refgenie)

    def lookup(parsed):
        if cmd.genome_digest:
            return refgenie.servers.seek(
                genome_digest=genome_ref(None, cmd.genome_digest),
                asset_group_name=parsed.asset_group,
                asset_name=parsed.asset,
                seek_key=parsed.seek_key,
                server_urls=server_urls,
            )
        return refgenie.servers.seek_components(parsed, server_urls=server_urls)

    run_over_paths(
        cmd.asset_registry_paths,
        lambda raw: _lookup_one(refgenie, raw, lookup, cmd.genome_digest),
    )


def _print_genome_id(cmd, refgenie, genome_digest) -> None:
    """Print a local genome's digest, validating its store or adding metadata if asked."""
    from refgenie.exceptions import RefgenieError

    if cmd.validate_store:
        # The store backend is pluggable (local filesystem, HTTP, S3), so the
        # transport errors it can raise are open-ended. Catch broadly, but
        # always terminate non-zero -- never fall through to a success exit.
        try:
            metadata = refgenie.store_router.get_collection_level2(genome_digest)
        except Exception as e:
            fail(f"Store validation failed for {genome_digest}: {e}")
        if metadata is None:
            fail(f"Genome {genome_digest} has no sequences in store", EXIT_NOT_FOUND)

    if cmd.verbose:
        try:
            metadata = refgenie.genome.get_metadata(genome_digest)
            rprint(f"digest: {metadata['digest']}")
            rprint(f"sequences: {metadata['n_sequences']}")
            rprint(f"total_length: {metadata['total_length']:,}")
            rprint(f"source: {metadata['source']}")
        except (RefgenieError, KeyError) as e:
            # Soft fallback: the digest itself is the command's output and
            # is already known; only the extra metadata lines are lost.
            logger.warning(f"Could not get metadata: {e}")
            rprint(genome_digest)
    else:
        rprint(genome_digest)


def _id_by_digest(cmd, refgenie) -> None:
    """``id --genome-digest D``: the genome locally, else on remote servers with --remote."""
    genome_digest = genome_ref(None, cmd.genome_digest)
    if refgenie.genome.exists(genome_digest):
        _print_genome_id(cmd, refgenie, genome_digest)
    elif cmd.remote:
        if not refgenie.servers.find_collection(genome_digest):
            fail(
                f"Digest '{genome_digest}' not found on any subscribed remote server.",
                EXIT_NOT_FOUND,
            )
        rprint(genome_digest)
    else:
        fail(
            f"Genome digest '{genome_digest}' is not known here. "
            f"Add --remote to look it up on subscribed servers.",
            EXIT_NOT_FOUND,
        )


def handle_id(cmd, refgenie) -> None:
    from refgenie.exceptions import MissingAliasError, RefgenieError
    from refgenie.models import GenomeAlias, GenomeDigest

    # Handle --info flag (digest-to-info lookup)
    if cmd.info:
        for digest_str in cmd.asset_registry_paths:
            try:
                digest = GenomeDigest(digest_str)
                metadata = refgenie.genome.get_metadata(digest)
                aliases = refgenie.alias.get_for_genome(genome_digest=digest)
                rprint(f"digest: {metadata['digest']}")
                rprint(f"aliases: {', '.join(aliases) if aliases else '(none)'}")
                rprint(f"sequences: {metadata['n_sequences']}")
                rprint(f"total_length: {metadata['total_length']:,}")
                rprint(f"source: {metadata['source']}")
                if metadata.get("remote_url"):
                    rprint(f"remote_url: {metadata['remote_url']}")
            except (RefgenieError, KeyError, ValueError) as e:
                fail(f"Cannot get info for digest '{digest_str}': {e}", EXIT_NOT_FOUND)
        return

    if not cmd.asset_registry_paths:
        if cmd.genome_digest:
            _id_by_digest(cmd, refgenie)
            return
        fail("Give a genome alias, an asset path, or --genome-digest.", EXIT_INVALID_INPUT)

    for path_str in cmd.asset_registry_paths:
        parsed = refgenie.parse_asset_registry_path(path_str)

        if parsed.genome is None and not cmd.genome_digest:
            # Genome-only input: a positional name is always an alias.
            name = parsed.asset_group
            try:
                genome_digest = refgenie.alias.resolve(GenomeAlias(name))
            except MissingAliasError:
                if cmd.remote:
                    fail(
                        f"'{name}' is not a known local alias. Remote lookup takes a digest: "
                        f"refgenie id --genome-digest <digest> --remote",
                        EXIT_INVALID_INPUT,
                    )
                fail(
                    f"'{name}' is not a known genome. "
                    f"For asset digests, use: refgenie id <genome>/{name}",
                    EXIT_NOT_FOUND,
                )
            _print_genome_id(cmd, refgenie, genome_digest)

        else:
            # Asset path - return asset digest. The genome is the path's alias,
            # or --genome-digest.
            genome_digest = genome_digest_from(refgenie, parsed.genome, cmd.genome_digest)
            asset_name = parsed.asset or refgenie.asset.group.get_default(
                genome_digest=genome_digest,
                asset_group_name=parsed.asset_group,
            )
            rprint(
                refgenie.asset.get(
                    genome_digest=genome_digest,
                    asset_group_name=parsed.asset_group,
                    asset_name=asset_name,
                ).digest
            )


def handle_compare(cmd, refgenie) -> None:
    count = len(cmd.genomes) + len(cmd.genome_digest or [])
    if count != 2:
        fail(
            f"compare takes exactly two genomes (aliases, or digests with --genome-digest); "
            f"got {count}.",
            EXIT_INVALID_INPUT,
        )
    digests = [genome_digest_from(refgenie, alias, None) for alias in cmd.genomes] + [
        genome_digest_from(refgenie, None, digest) for digest in cmd.genome_digest or []
    ]
    rprint(refgenie.genome.compare(*digests))
