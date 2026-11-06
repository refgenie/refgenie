"""Filesystem helpers for the alias symlink trees (no database access)."""

import os
import shutil
from inspect import signature
from pathlib import Path
from collections.abc import Callable

from refgenie.const import ALIAS_DIR
from refgenie.logger import logger


def _replace_genome_digest(text: str, genome_digest: str, alias_name: str) -> str:
    """Replace the genome digest with the alias name in a string."""
    return text.replace(genome_digest, alias_name)


def get_symlink_paths(
    alias_folder: Path,
    aliases: list[str],
    asset_group_name: str | None = None,
    asset_name: str | None = None,
) -> dict[str, Path]:
    """
    Get paths to the alias directories for the given aliases and optional asset.

    Paths go as deep as the arguments allow: genome, group, or asset level.
    Group cleanup uses the group level, since a default asset name may be stale.

    Args:
        alias_folder: The base alias folder path.
        aliases: List of alias names to create paths for.
        asset_group_name: The name of the asset group (optional).
        asset_name: The name of the asset (optional).

    Returns:
        dict[str, Path]: Mapping of alias names to their directory paths.
    """
    if len(aliases) == 0:
        return {}
    if asset_group_name:
        if asset_name is None:
            return {alias: alias_folder / alias / asset_group_name for alias in aliases}
        return {alias: alias_folder / alias / asset_group_name / asset_name for alias in aliases}
    return {alias: alias_folder / alias for alias in aliases}


def get_build_paths(
    genome_folder: Path,
    aliases: list[str],
    asset_group_name: str | None = None,
    asset_name: str | None = None,
) -> dict[str, Path]:
    """
    Get paths to the build directories for the given aliases and optional asset.

    Like :func:`get_symlink_paths`, but rooted at the ``builds/`` tree.

    Args:
        genome_folder: The refgenie genome folder.
        aliases: List of alias names to create paths for.
        asset_group_name: The name of the asset group (optional).
        asset_name: The name of the asset (optional).

    Returns:
        dict[str, Path]: Mapping of alias names to their build directory paths.
    """
    from refgenie.utils.paths import get_build_dir

    if len(aliases) == 0:
        return {}
    return {
        alias: get_build_dir(
            genome_folder=genome_folder,
            genome_name=alias,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        for alias in aliases
    }


def alias_owned_paths(genome_folder: Path, alias_name: str) -> list[Path]:
    """
    The on-disk trees owned by a genome alias.

    ``alias/<name>/`` and ``builds/<name>/``. Both alias backends use this
    list and queue the removal until the alias delete has committed.

    Args:
        genome_folder: The refgenie genome folder.
        alias_name: The alias whose trees to enumerate.

    Returns:
        list[Path]: The alias directory and the build directory, in that order.
    """
    from refgenie.utils.paths import get_build_dir

    return [
        genome_folder / ALIAS_DIR / alias_name,
        get_build_dir(genome_folder=genome_folder, genome_name=alias_name),
    ]


def remove_alias_files(genome_folder: Path, alias_name: str) -> None:
    """
    Remove the trees :func:`alias_owned_paths` names, immediately.

    Only for callers with no transaction to order against.

    Args:
        genome_folder: The refgenie genome folder.
        alias_name: The alias whose files should be removed.
    """
    for path in alias_owned_paths(genome_folder=genome_folder, alias_name=alias_name):
        logger.info(f"Deleting alias-owned files: {path}")
        shutil.rmtree(path, ignore_errors=True)


def create_alias_symlinks(
    src_path: Path,
    target_paths_mapping: dict[str, Path],
    genome_digest: str,
    link_fun: Callable[[str, str], None] = lambda t, s: os.symlink(t, s),  # noqa: E731 - lambda documents (target, symlink) arg order
) -> list[Path]:
    """
    Create symlinks from source path to alias directories.

    Recreates the source tree in each target path with symlinks to the files,
    renaming the genome digest in file and directory names to that target's
    alias. Idempotent: existing symlinks are repointed; a real file raises.

    Args:
        src_path: The source directory to symlink from.
        target_paths_mapping: Mapping of alias names to their target directory paths.
        genome_digest: The genome digest (used for string replacement in filenames).
        link_fun: Function to use for linking (default: os.symlink).
            Must accept two arguments: (target, destination).

    Returns:
        list[Path]: List of created alias directory paths.

    Raises:
        TypeError: If link_fun is not callable.
        FileExistsError: If a destination path exists and is not a symlink.
    """

    def _rpl(x: str, alias: str) -> str:
        return _replace_genome_digest(x, genome_digest, alias)

    if not callable(link_fun):
        raise TypeError("link_fun must be callable")
    try:
        params = signature(link_fun).parameters
        required = [p for p in params.values() if p.default is p.empty]
        if len(required) != 2:
            raise TypeError(
                "link_fun must accept exactly 2 required arguments (target, destination)"
            )
    except (ValueError, TypeError):
        pass  # Can't inspect signature (e.g., builtins) - let it fail at call time if wrong

    created = []
    for alias, path in target_paths_mapping.items():
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        for root, dirs, files in os.walk(src_path):
            appendix = os.path.relpath(root, src_path)
            for directory in dirs:
                os.makedirs(os.path.join(path, appendix, _rpl(directory, alias)), exist_ok=True)
            for file in files:
                target = os.path.join(root, file)  # The actual file path to link to
                new_path = os.path.join(path, appendix, _rpl(file, alias))  # The link path
                if os.path.islink(new_path):
                    os.unlink(new_path)
                elif os.path.exists(new_path):
                    raise FileExistsError(
                        f"Refusing to replace non-symlink in alias tree: {new_path}"
                    )
                link_fun(target, new_path)
        created.append(path)

    if created:
        logger.info(f"Created alias directories: {', '.join(p.as_posix() for p in created)}")

    return created
