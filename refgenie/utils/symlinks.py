"""
Symlink utility functions for creating alias directory structures.

These are pure filesystem utilities with no database access.
"""

import os
import shutil
from inspect import signature
from pathlib import Path
from collections.abc import Callable

from refgenie.const import ALIAS_DIR
from refgenie.logger import logger


def _replace_genome_digest(text: str, genome_digest: str, alias_name: str) -> str:
    """
    Replace genome digest with human-readable genome alias in a string.

    Args:
        text: String to process.
        genome_digest: The genome digest to replace.
        alias_name: The alias name to replace the digest with.

    Returns:
        str: Processed string with digest replaced by alias.
    """
    return text.replace(genome_digest, alias_name)


def get_symlink_paths(
    alias_folder: Path,
    aliases: list[str],
    asset_group_name: str | None = None,
    asset_name: str | None = None,
) -> dict[str, Path]:
    """
    Get paths to the alias directories for the given aliases and optional asset.

    Paths resolve as deeply as the given arguments allow: genome-level with no
    asset group, group-level with a group but no asset name, asset-level with
    both. The group-level form is what group cleanup needs — it must not
    substitute a default asset name, which at removal time is stale.

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

    Mirrors :func:`get_symlink_paths`, but rooted at the ``builds/`` tree rather
    than the alias tree. Used by the removal paths to tear down build
    bookkeeping alongside the alias directories.

    Args:
        genome_folder: The refgenie genome folder.
        aliases: List of alias names to create paths for.
        asset_group_name: The name of the asset group (optional).
        asset_name: The name of the asset (optional).

    Returns:
        dict[str, Path]: Mapping of alias names to their build directory paths.
    """
    from refgenie.utils.build import get_build_dir

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

    An alias owns exactly two trees: ``alias/<name>/``, the name-addressed view
    of the data directory, and ``builds/<name>/``, its build bookkeeping.
    Neither outlives the alias that names them.

    Both alias backends use it — the SQL-backed ``AliasManager`` via the
    ``Alias`` before_delete handler, and the store-backed ``StoreAliasManager``
    directly. Adding a third alias-keyed tree means editing this function, not
    hunting for every backend that deletes an alias.

    It returns paths rather than removing them: removal is deferred until the
    catalog that owns the alias has committed its removal, so the two backends
    queue these paths rather than acting on them.

    Args:
        genome_folder: The refgenie genome folder.
        alias_name: The alias whose trees to enumerate.

    Returns:
        list[Path]: The alias directory and the build directory, in that order.
    """
    from refgenie.utils.build import get_build_dir

    return [
        genome_folder / ALIAS_DIR / alias_name,
        get_build_dir(genome_folder=genome_folder, genome_name=alias_name),
    ]


def remove_alias_files(genome_folder: Path, alias_name: str) -> None:
    """
    Remove the trees :func:`alias_owned_paths` names, immediately.

    Only for callers with no transaction to order against — the alias backends
    defer instead.

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

    Walks the source directory tree and recreates the structure in each target path,
    creating symbolic links to the actual files instead of copying them.

    ``target_paths_mapping`` is the sole source of alias names: each target
    directory has its file and directory names rewritten with the alias that
    owns it, so a genome with several aliases gets each tree named consistently
    rather than all of them named for one arbitrary alias.

    Creation is idempotent. A destination that is already a symlink is
    repointed; a destination that is a real file is left alone and raises,
    since that is data in the alias tree rather than a link this function owns.

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
    # TODO: consider rewriting this function to use pathlib once Path.walk() is available (3.12+)

    def _rpl(x: str, alias: str) -> str:
        """Replace genome digest with alias name in string."""
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
