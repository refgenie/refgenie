"""Field descriptions and handler helpers shared across command families."""

import sys

from refgenie.cli.errors import (
    EXIT_GENERAL_ERROR,
    EXIT_INVALID_INPUT,
    EXIT_NOT_FOUND,
    EXIT_OK,
    fail,
)
from refgenie.logger import logger

# --- Shared field descriptions ---

_GENOME_DESCRIPTION = "Genome alias, e.g. mm10."
_GENOME_DIGEST_DESCRIPTION = (
    "Genome digest, used in place of an alias. Reaches a genome with no alias. "
    "With a registry path, leave the genome out of the path (e.g. fasta, not hg38/fasta)."
)
_ASSET_REGISTRY_PATHS_DESCRIPTION = (
    "One or more registry path strings that identify assets "
    "(e.g. hg38/fasta or hg38/fasta:default)."
)
_ASSET_REGISTRY_PATHS_SEEK_DESCRIPTION = (
    "One or more registry path strings that identify assets "
    "(e.g. hg38/fasta or hg38/fasta:default or hg38/fasta.fai:default)."
)


# Inlined on each command rather than a mixin, which would reorder flags in --help.
_GENOME_SERVER_TRANSIENT_DESCRIPTION = (
    "One or more URLs to use. This information will not persist in the genome config file."
)
_APPEND_SERVER_DESCRIPTION = "Whether the provided servers should be appended to the list."


# --- Shared handler helpers ---


def _validate_registry_paths(parsed_paths, raw_paths, genome_digest: str | None = None):
    """Exit with a friendly error if any parsed registry path is missing its parts.

    A path needs a genome unless ``--genome-digest`` names it, and must not
    have one when it does.
    """
    for i, parsed in enumerate(parsed_paths):
        if not parsed.asset_group or not (parsed.genome or genome_digest):
            fail(
                f"Invalid asset path: '{raw_paths[i]}'. "
                "Expected format: genome/asset_group[.seek_key][:asset] (e.g. rCRSd/fasta)",
                EXIT_INVALID_INPUT,
            )
        genome_ref(parsed.genome, genome_digest)


def genome_ref(path_genome: str | None, genome_digest: str | None):
    """The genome a command names, typed by how it was named.

    The genome in a registry path (or a ``-g`` value) is an alias;
    ``--genome-digest`` is a digest. Nothing here guesses which one a string is.

    Args:
        path_genome: The alias: the genome part of a registry path, if any.
        genome_digest: The ``--genome-digest`` value, if any.

    Returns:
        GenomeAlias | GenomeDigest: The typed identifier.
    """
    from refgenie.models import GenomeAlias, GenomeDigest

    if path_genome and genome_digest:
        fail(
            f"Name the genome once: by its alias ('{path_genome}') or with --genome-digest, "
            f"not both.",
            EXIT_INVALID_INPUT,
        )
    try:
        if genome_digest:
            return GenomeDigest(genome_digest)
        if path_genome:
            return GenomeAlias(path_genome)
    except ValueError as e:
        fail(str(e), EXIT_INVALID_INPUT)
    fail(
        "No genome given. Put its alias in the registry path (e.g. hg38/fasta) "
        "or pass --genome-digest.",
        EXIT_INVALID_INPUT,
    )


def genome_digest_from(refgenie, path_genome: str | None, genome_digest: str | None):
    """The digest of the genome a command names: the flag as is, or the alias resolved.

    Raises:
        MissingAliasError: If the alias is not known here.
    """
    from refgenie.models import GenomeDigest

    ref = genome_ref(path_genome, genome_digest)
    return ref if isinstance(ref, GenomeDigest) else refgenie.alias.resolve(ref)


def check_push_remotes(refgenie, refs: list[str] | None) -> None:
    """Exit unless every ``--push-to`` remote exists, before any build or stage starts."""
    from refgenie.exceptions import MissingRemoteError

    for ref in refs or []:
        try:
            refgenie.remote.get(ref)
        except MissingRemoteError as e:
            fail(f"--push-to: {e}", EXIT_NOT_FOUND)


def genome_digests_from(
    refgenie, aliases: list[str] | None, genome_digests: list[str] | None
) -> list | None:
    """Digests for a genome filter given as ``-g`` aliases and ``--genome-digest`` digests.

    Returns None when neither is given, meaning "no filter".

    Raises:
        MissingAliasError: If an alias is not known here.
    """
    if not aliases and not genome_digests:
        return None
    return [genome_digest_from(refgenie, a, None) for a in aliases or []] + [
        genome_digest_from(refgenie, None, d) for d in genome_digests or []
    ]


def data_channel_hint(retry_command: str | None = None) -> str:
    """Return the standard guidance for registering and syncing a data channel.

    Fresh installs have no recipes or asset classes until a channel is synced.

    Args:
        retry_command: If given, appended as a final line (e.g.
            "refgenie build hg38/fasta").
    """
    lines = [
        "Register and sync a data channel first, then retry:",
        "  refgenie data-channel add refgenie https "
        "https://refgenie.github.io/refgenie-registry/index.yaml",
        "  refgenie data-channel sync refgenie --exists-ok",
    ]
    if retry_command:
        lines.append(f"  {retry_command}")
    return "\n".join(lines)


def no_data_channel_synced(refgenie) -> bool:
    """True when no recipes or asset classes are registered at all.

    Used to decide whether a missing recipe or asset class means "sync a data
    channel" rather than "that name does not exist".
    """
    return not any(refgenie.recipe.list_all()) and not any(refgenie.asset_class.list_all())


def resolve_transient_servers(cmd, refgenie) -> list[str] | None:
    """Resolve --genome-server/--append-server into a per-invocation URL list.

    Returns None when no --genome-server was given, meaning "use subscriptions".
    Never writes to the database.
    """
    if not cmd.genome_server:
        return None
    urls = list(cmd.genome_server)
    if cmd.append_server:
        urls = refgenie.servers.subscriptions() + urls
    return list(dict.fromkeys(urls))


def _lookup_one(refgenie, raw: str, lookup, genome_digest: str | None = None) -> int:
    """Resolve one registry path with ``lookup`` and print the result.

    With ``genome_digest`` (the ``--genome-digest`` flag), the path names no
    genome; the flag does.

    Returns an exit code: EXIT_INVALID_INPUT, EXIT_NOT_FOUND, EXIT_GENERAL_ERROR
    or EXIT_OK.
    """
    from refgenie.exceptions import (
        MissingAliasError,
        MissingAssetError,
        MissingAssetGroupError,
        MissingGenomeError,
        MissingSeekKeyError,
        RefgenieError,
    )

    parsed = refgenie.parse_asset_registry_path(raw)
    if not parsed.asset_group or not (parsed.genome or genome_digest):
        logger.error(
            f"Invalid asset path: '{raw}'. "
            "Expected format: genome/asset_group[.seek_key][:asset] (e.g. rCRSd/fasta)"
        )
        return EXIT_INVALID_INPUT
    genome_ref(parsed.genome, genome_digest)
    try:
        sys.stdout.write(str(lookup(parsed)) + "\n")
    except (
        MissingAliasError,
        MissingGenomeError,
        MissingAssetGroupError,
        MissingAssetError,
        MissingSeekKeyError,
    ) as e:
        logger.error(f"{raw}: {e}")
        return EXIT_NOT_FOUND
    except (ValueError, FileNotFoundError, RefgenieError) as e:
        logger.error(f"{raw}: {e}")
        return EXIT_GENERAL_ERROR
    return EXIT_OK
