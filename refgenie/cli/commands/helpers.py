"""Pieces shared across command families: the field descriptions that recur in
their models, and the helpers their handlers call.

Separate from ``framework.py``, which holds parser machinery only.
"""

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


# Transient server-selection fields recur on seekr, listr, and populater, but
# are inlined (not a mixin): pydantic orders base-class fields before subclass
# fields, and a mixin would reorder each command's flags in --help.
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
                f"Expected format: genome/asset_class (e.g. rCRSd/fasta)",
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

    Fresh installs ship with no recipes or asset classes -- these only arrive
    after a data channel is registered and synced. This is the single source
    of truth for the hint shown when a command fails for that reason; `pull`,
    `build`, and `genome init` all call this instead of hand-rolling the
    message so the guidance stays byte-identical everywhere it appears.

    Pure string builder: no logging, no `sys.exit`. Callers control their own
    control flow and decide whether/how to emit the returned text.

    Args:
        retry_command: If given, appended as a final line so the user can
            re-run the command that originally failed (e.g.
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

    A `MissingRecipeError`/`MissingAssetClassError` raised while both
    registries are completely empty is the signature of an unsynced data
    channel (show the hint). The same error when other recipes/classes exist
    means the user asked for one that genuinely is not in their synced
    channel -- do not misdirect them to sync a channel in that case.
    """
    return not any(refgenie.recipe.list_all()) and not any(refgenie.asset_class.list_all())


def resolve_transient_servers(cmd, refgenie) -> list[str] | None:
    """Resolve --genome-server/--append-server into a per-invocation URL list.

    Returns None when no --genome-server was given, meaning "use subscriptions".
    Never writes to the database: these flags are documented as transient.
    """
    if not cmd.genome_server:
        return None
    urls = list(cmd.genome_server)
    if cmd.append_server:
        urls = refgenie.servers.subscriptions() + urls
    # De-duplicate, preserving order.
    return list(dict.fromkeys(urls))


def _lookup_one(refgenie, raw: str, lookup, genome_digest: str | None = None) -> int:
    """Resolve one registry path with ``lookup`` and print the result.

    With ``genome_digest`` (the ``--genome-digest`` flag), the path names no
    genome; the flag does.

    Returns an exit code: EXIT_INVALID_INPUT for a malformed path,
    EXIT_NOT_FOUND for a missing genome/asset/seek key, EXIT_GENERAL_ERROR for
    anything else, EXIT_OK on success.
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
            f"Invalid asset path: '{raw}'. Expected format: genome/asset_class (e.g. rCRSd/fasta)"
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
