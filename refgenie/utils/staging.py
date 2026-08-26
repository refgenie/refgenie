"""Staged asset path conventions.

The staged archive path is content-addressed by asset digest:
    {genome_stage_folder}/{genome_digest}/{group_name}/{asset_digest}.tgz

This module centralizes that convention. All code that constructs staged
archive paths should use these helpers to ensure consistency.
"""

from pathlib import Path


def staged_archive_relpath(genome_digest: str, group_name: str, asset_digest: str) -> str:
    """Return the relative path for a staged archive.

    Args:
        genome_digest: The genome's digest (e.g. "abc123...").
        group_name: The asset group name (e.g. "bowtie2_index").
        asset_digest: The asset's content digest.

    Returns:
        Relative path string: "{genome_digest}/{group_name}/{asset_digest}.tgz"
    """
    return f"{genome_digest}/{group_name}/{asset_digest}.tgz"


def staged_archive_path(
    genome_stage_folder: Path | str,
    genome_digest: str,
    group_name: str,
    asset_digest: str,
) -> Path:
    """Return the full path for a staged archive.

    Args:
        genome_stage_folder: Base staging folder.
        genome_digest: The genome's digest.
        group_name: The asset group name.
        asset_digest: The asset's content digest.

    Returns:
        Absolute path: {genome_stage_folder}/{genome_digest}/{group_name}/{asset_digest}.tgz
    """
    return Path(genome_stage_folder) / staged_archive_relpath(genome_digest, group_name, asset_digest)
