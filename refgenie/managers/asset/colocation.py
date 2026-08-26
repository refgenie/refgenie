"""
Colocation utilities for managing parent-asset symlinks in child asset directories.

Recipes can declare `colocate` requirements on input assets, specifying which
parent files should be symlinked into the child's output directory. This module
handles creating those symlinks (at build time and after pull) and identifying
them for archive exclusion.
"""

import os
from pathlib import Path

from refgenie.db.tables import Asset
from refgenie.logger import logger


def link_colocated_file(
    output_folder: Path,
    genome_folder: Path,
    parent_asset: Asset,
    source_key: str,
    dest: str | None = None,
) -> str | None:
    """Create one relative symlink to a parent asset's file.

    The single place the build path and the pull path agree on: resolve the
    parent's seek key to a real file, name the link, and point it at the parent
    relatively so the whole genome folder stays movable.

    Args:
        output_folder: The directory the link is created in.
        genome_folder: The root genome folder, which seek key paths are relative to.
        parent_asset: The asset holding the file being linked to.
        source_key: The parent's seek key naming that file.
        dest: The link name. Defaults to the parent file's own name.

    Returns:
        The created link's filename, or None if nothing was created (the parent
        file is missing, the seek key is unknown, or the name is already taken).
    """
    rel_path = parent_asset.seek_keys_dict.get(source_key)
    if rel_path is None:
        logger.warning(
            f"Seek key '{source_key}' not found in parent asset. "
            f"Available: {list(parent_asset.seek_keys_dict.keys())}"
        )
        return None

    parent_file = genome_folder / rel_path
    if not parent_file.exists():
        logger.warning(f"Parent file does not exist: {parent_file}")
        return None

    dest_name = dest or parent_file.name
    link_path = output_folder / dest_name
    if link_path.exists() or link_path.is_symlink():
        logger.debug(f"Colocation target already exists: {link_path}")
        return None

    rel_target = os.path.relpath(parent_file, output_folder)
    link_path.symlink_to(rel_target)
    logger.info(f"Created colocation symlink: {link_path} -> {rel_target}")
    return dest_name


def recreate_colocation_symlinks(
    output_folder: Path,
    genome_folder: Path,
    colocate_metadata: list[dict[str, str]],
    parent_assets: dict[str, Asset | None],
) -> list[str]:
    """Recreate colocation symlinks after a pull extraction.

    The build path works from the recipe's ``input_assets``; a pulled asset has
    no recipe on hand, only the flat ``colocate`` metadata stored on it. Both
    end at the same symlink, made by `link_colocated_file`.

    Args:
        output_folder: The extracted asset directory.
        genome_folder: The root genome folder.
        colocate_metadata: Entries with ``parent_asset_group``, ``source_key``,
            and optional ``dest``.
        parent_assets: Parent asset group name -> the resolved parent asset.

    Returns:
        List of filenames created in output_folder.
    """
    created_files = []
    for entry in colocate_metadata:
        parent_group = entry.get("parent_asset_group")
        source_key = entry.get("source_key")
        if not parent_group or not source_key:
            continue
        parent_asset = parent_assets.get(parent_group)
        if parent_asset is None:
            logger.warning(f"Parent asset '{parent_group}' not resolved, skipping colocation")
            continue
        created = link_colocated_file(
            output_folder, genome_folder, parent_asset, source_key, entry.get("dest")
        )
        if created is not None:
            created_files.append(created)
    return created_files


def create_colocation_symlinks(
    output_folder: Path,
    genome_folder: Path,
    input_assets: dict | None,
    resolved_assets: dict[str, Asset | None] | None,
) -> list[str]:
    """Create relative symlinks for colocated parent files.

    Walks the recipe's input_assets looking for `colocate` declarations, resolves
    each parent's source_key to an absolute path, and creates a relative symlink
    in output_folder.

    Args:
        output_folder: The child asset's output directory.
        genome_folder: The root genome folder (to resolve relative seek key paths).
        input_assets: The recipe's input_assets dict, e.g.
            {"fasta": {"asset_class": "fasta", "default": "fasta",
                       "colocate": [{"source_key": "fasta", "dest": "genome.fa"}]}}
        resolved_assets: Mapping of input asset names to resolved Asset objects.

    Returns:
        List of filenames created in output_folder (for archive exclusion).
    """
    if not input_assets or not resolved_assets:
        return []

    created_files = []
    output_folder.mkdir(parents=True, exist_ok=True)

    for input_name, input_spec in input_assets.items():
        colocate_list = input_spec.get("colocate")
        if not colocate_list:
            continue

        parent_asset = resolved_assets.get(input_name)
        if parent_asset is None:
            logger.warning(f"Parent asset '{input_name}' not resolved, skipping colocation")
            continue

        for entry in colocate_list:
            source_key = entry.get("source_key")
            if not source_key:
                logger.warning(f"No source_key in colocate entry for '{input_name}'")
                continue
            created = link_colocated_file(
                output_folder, genome_folder, parent_asset, source_key, entry.get("dest")
            )
            if created is not None:
                created_files.append(created)

    return created_files


def get_colocation_filenames(
    input_assets: dict | None,
    resolved_assets: dict[str, Asset | None] | None,
) -> list[str]:
    """Return the list of filenames that colocation would create.

    Used to determine which files to exclude from archives without
    needing to inspect the filesystem.

    Args:
        input_assets: The recipe's input_assets dict.
        resolved_assets: Mapping of input asset names to resolved Asset objects.

    Returns:
        List of filenames that would be created by create_colocation_symlinks.
    """
    if not input_assets or not resolved_assets:
        return []

    filenames = []
    for input_name, input_spec in input_assets.items():
        colocate_list = input_spec.get("colocate")
        if not colocate_list:
            continue

        parent_asset = resolved_assets.get(input_name)
        if parent_asset is None:
            continue

        parent_seek_keys = parent_asset.seek_keys_dict

        for entry in colocate_list:
            source_key = entry.get("source_key")
            if not source_key:
                continue

            if dest := entry.get("dest"):
                filenames.append(dest)
            elif (rel_path := parent_seek_keys.get(source_key)) is not None:
                filenames.append(Path(rel_path).name)

    return filenames


def get_colocation_metadata(
    input_assets: dict | None,
) -> list[dict[str, str]] | None:
    """Extract colocation metadata from recipe input_assets for storage on the Asset.

    Returns a list of dicts suitable for JSON storage, containing all info
    needed to recreate symlinks after pull.

    Args:
        input_assets: The recipe's input_assets dict.

    Returns:
        List of colocation metadata dicts, or None if no colocation is declared.
    """
    if not input_assets:
        return None

    metadata = []
    for input_name, input_spec in input_assets.items():
        colocate_list = input_spec.get("colocate")
        if not colocate_list:
            continue

        default_asset_group = input_spec.get("default", input_name)
        for entry in colocate_list:
            source_key = entry.get("source_key")
            if not source_key:
                continue
            item = {
                "parent_asset_group": default_asset_group,
                "source_key": source_key,
            }
            if dest := entry.get("dest"):
                item["dest"] = dest
            metadata.append(item)

    return metadata or None
