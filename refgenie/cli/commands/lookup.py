"""The commands that resolve one thing and print it -- `seek`, `seekr`, `id`,
`compare`: models and handlers.
"""

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg
from rich import print as rprint

from refgenie.cli.commands.framework import CliList
from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _ASSET_REGISTRY_PATHS_SEEK_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    _lookup_one,
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


class IdModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description="One or more registry paths: genome name (e.g. hg38) for "
        "genome digest, or asset path (e.g. hg38/fasta or hg38/fasta:default) "
        "for asset digest.",
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
        description="Query subscribed seqcolapi servers if genome is not found locally.",
        validation_alias=AliasChoices("remote"),
    )
    info: bool = Field(
        default=False,
        description="Given a digest, show its aliases and metadata.",
        validation_alias=AliasChoices("info"),
    )


class CompareModel(BaseModel):
    genome1: CliPositionalArg[str] = Field(
        description="First genome for compatibility check.",
    )
    genome2: CliPositionalArg[str] = Field(
        description="Second genome for compatibility check.",
    )


def handle_seek(cmd, refgenie) -> None:
    def lookup(parsed):
        return refgenie.asset._seek_by_components(
            parsed,
            force_exists=getattr(cmd, "check_exists", False),
            abs_path=getattr(cmd, "abs", False),
        )

    run_over_paths(
        cmd.asset_registry_paths,
        lambda raw: _lookup_one(refgenie, raw, lookup),
    )


def handle_seekr(cmd, refgenie) -> None:
    server_urls = resolve_transient_servers(cmd, refgenie)

    def lookup(parsed):
        return refgenie.asset.seek_remote(
            genome_name=parsed.genome,
            asset_group_name=parsed.asset_group,
            asset_name=parsed.asset,
            seek_key=parsed.seek_key,
            server_urls=server_urls,
        )

    run_over_paths(
        cmd.asset_registry_paths,
        lambda raw: _lookup_one(refgenie, raw, lookup),
    )


def handle_id(cmd, refgenie) -> None:
    import re

    from refgenie.exceptions import MissingAliasError, RefgenieError

    # Handle --info flag (digest-to-info lookup)
    if cmd.info:
        for digest_str in cmd.asset_registry_paths:
            try:
                metadata = refgenie.genome.get_metadata(digest_str)
                aliases = refgenie.alias.get_for_genome(genome_digest=digest_str)
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

    for path_str in cmd.asset_registry_paths:
        parsed = refgenie.parse_asset_registry_path(path_str)

        if parsed.genome is None:
            # Genome-only input
            name = parsed.asset_group

            # Try local resolution first
            try:
                genome_digest = refgenie.alias.resolve(name)
            except MissingAliasError:
                genome_digest = None

            # Handle --remote flag
            # Remote lookup requires a digest, not a name.
            # Seqcolapi collections do not have name fields matching refgenie aliases.
            if genome_digest is None and cmd.remote:
                if re.fullmatch(r"[a-zA-Z0-9_-]{32,}", name):
                    remote_info = refgenie.check_remote_digest(name)
                    if remote_info:
                        rprint(name)
                        return
                    else:
                        fail(
                            f"Digest '{name}' not found on any subscribed remote server.",
                            EXIT_NOT_FOUND,
                        )
                else:
                    fail(
                        f"'{name}' is not a known local alias and remote lookup requires a digest, "
                        f"not a name. Resolve the alias locally first with 'refgenie id {name}', "
                        f"then pass the resulting digest with --remote.",
                        EXIT_INVALID_INPUT,
                    )

            if genome_digest is None:
                fail(
                    f"'{name}' is not a known genome. "
                    f"For asset digests, use: refgenie id <genome>/{name}",
                    EXIT_NOT_FOUND,
                )

            # Validate store if requested
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

            # Output
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

        else:
            # Full asset path - return asset digest
            asset_name = parsed.asset or refgenie.asset.get_default(
                genome_name=parsed.genome,
                asset_group_name=parsed.asset_group,
            )
            rprint(
                refgenie.asset.get(
                    genome_name=parsed.genome,
                    asset_group_name=parsed.asset_group,
                    asset_name=asset_name,
                ).digest
            )


def handle_compare(cmd, refgenie) -> None:
    digest_a = refgenie.alias.resolve(cmd.genome1)
    digest_b = refgenie.alias.resolve(cmd.genome2)
    rprint(refgenie.genome.compare(digest_a, digest_b))
