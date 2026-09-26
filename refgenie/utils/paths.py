"""The ``builds/`` bookkeeping-directory formula, shared by build, symlinks, and job logs."""

from pathlib import Path

from refgenie.const import BUILDS_DIR


def get_build_dir(
    genome_folder: Path,
    genome_name: str,
    asset_group_name: str | None = None,
    asset_name: str | None = None,
) -> Path:
    """
    Get the build bookkeeping directory for a build invocation.

    This is the single formula for the ``builds/`` layout; every call site
    derives its path from here rather than joining the components itself.
    The tree is keyed by build invocation (alias/genome name, group, asset),
    NOT by any asset addressing scheme — it must be computable before the
    build runs, and it must not follow asset content when content moves.

    Args:
        genome_folder: The refgenie genome folder.
        genome_name: The genome alias name (or a ``{genome_name}`` template
            placeholder for snakemake substitution).
        asset_group_name: The name of the asset group (optional).
        asset_name: The name of the asset (optional).

    Returns:
        Path: The build directory, resolved as deeply as the given arguments allow.
    """
    build_dir = genome_folder / BUILDS_DIR / genome_name
    if asset_group_name is None:
        return build_dir
    build_dir = build_dir / asset_group_name
    if asset_name is None:
        return build_dir
    return build_dir / asset_name
