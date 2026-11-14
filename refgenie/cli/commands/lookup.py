"""The commands that resolve one thing and print it -- `seek`, `seekr`, `id`,
`compare`: models and handlers.
"""

import json
import sys

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg
from rich import print as rprint

from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _ASSET_REGISTRY_PATHS_SEEK_DESCRIPTION,
    _GENOME_DIGEST_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    _lookup_one,
    genome_digest_from,
    genome_ref,
    locate_genome,
    note,
    resolve_transient_servers,
)
from refgenie.cli.errors import EXIT_INVALID_INPUT, EXIT_NOT_FOUND, fail, run_over_paths
from refgenie.cli.framework import CliList
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
        description="Query subscribed refgenie servers if the --genome-digest genome "
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
        "exactly two genomes. A genome not registered here is compared on the "
        "subscribed server that has it.",
    )
    genome_digest: CliList | None = Field(
        default=None,
        description="Genome digests to compare, in place of aliases.",
        validation_alias=AliasChoices("genome-digest"),
    )
    as_json: bool = Field(
        default=False,
        description="Print the summary, the rename map and the GA4GH comparison as JSON.",
        validation_alias=AliasChoices("json"),
    )


def handle_seek(cmd, refgenie) -> None:
    force_exists = cmd.check_exists
    abs_path = cmd.abs

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
        lambda raw: _lookup_one(refgenie, raw, lookup, cmd.genome_digest, pull_hint=True),
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
    from refgenie.exceptions import MissingGenomeError, RefgenieError

    if cmd.validate_store:
        # Store backends (local, HTTP, S3) raise open-ended transport errors.
        try:
            metadata = refgenie.genome.level2(genome_digest)
        except MissingGenomeError as e:
            fail(f"Genome {genome_digest} has no sequences in store: {e}", EXIT_NOT_FOUND)
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
            # The digest is still printed; only the extra lines are lost.
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
        if not refgenie.servers.locate(genome_digest):
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
            # Asset path: print the asset digest.
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
    import httpx2

    from refgenie.utils.seqcol import compare_level2, summarize_comparison

    count = len(cmd.genomes) + len(cmd.genome_digest or [])
    if count != 2:
        fail(
            f"compare takes exactly two genomes (aliases, or digests with --genome-digest); "
            f"got {count}.",
            EXIT_INVALID_INPUT,
        )
    refs = [genome_ref(a, None) for a in cmd.genomes] + [
        genome_ref(None, d) for d in cmd.genome_digest or []
    ]
    (digest_a, server_a), (digest_b, server_b) = (locate_genome(refgenie, r) for r in refs)
    if server_a is None and server_b is None:
        side_a, side_b = refgenie.genome.compare_sides(digest_a, digest_b)
        comparison = compare_level2(side_a, side_b, digest_a, digest_b)
        summary, name_map = summarize_comparison(side_a, side_b, comparison, *map(str, refs))
    else:
        server_url = server_a or server_b
        try:
            comparison = refgenie.servers.client(server_url).genome.compare(digest_a, digest_b)
        except (httpx2.HTTPError, ValueError) as e:
            remote = [r for r, s in zip(refs, (server_a, server_b)) if s]
            fail(
                f"{server_url} could not compare them ({e}). Register "
                + ", ".join(f"'{r}'" for r in remote)
                + " here first (refgenie pull --init, metadata only), then compare again.",
                EXIT_NOT_FOUND,
            )
        summary, name_map = comparison.pop("summary"), comparison.pop("name_map")
    if cmd.as_json:
        payload = {"summary": summary, "name_map": name_map, "comparison": comparison}
        sys.stdout.write(json.dumps(payload, indent=2) + "\n")
        return
    lines = [summary] + [f"  {a} -> {b}" for a, b in name_map.items()]
    sys.stdout.write("\n".join(lines) + "\n")
    note("`refgenie compare --json` prints the full GA4GH comparison.")
