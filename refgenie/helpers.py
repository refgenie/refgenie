from __future__ import annotations

import errno
import os
from collections.abc import Iterable

from refgenconf import MissingRecipeError
from ubiquerg import is_writable

from .asset_build_packages import asset_build_packages
from .exceptions import MissingFolderError


def _parse_user_build_input(input: Iterable[Iterable[str]] | None) -> dict[str, str] | None:
    """Parse user input specification for build command.

    Used in build for specific parents and input parsing.

    Args:
        input: User command line input, formatted as
            [[fasta=txt, test=txt], ...].

    Returns:
        Mapping of input names to values, or None if input is empty.
    """
    lst = []
    for i in input or []:
        lst.extend(i)
    return (
        {x.split("=")[0]: x.split("=")[1] for x in lst if "=" in x}
        if lst is not None
        else lst
    )


def _single_folder_writeable(d: str) -> bool:
    return os.access(d, os.W_OK) and os.access(d, os.X_OK)


def _writeable(outdir: str | None, strict_exists: bool = False) -> bool:
    """Check whether an output directory is writeable.

    Args:
        outdir: Path to check. Defaults to current directory if None.
        strict_exists: Whether to raise an error if the path doesn't exist.

    Returns:
        True if the directory is writeable.

    Raises:
        MissingFolderError: If strict_exists is True and the path doesn't exist.
    """
    outdir = outdir or "."
    if os.path.exists(outdir):
        return _single_folder_writeable(outdir)
    elif strict_exists:
        raise MissingFolderError(outdir)
    return _writeable(os.path.dirname(outdir), strict_exists)


def _raise_missing_recipe_error(recipe: str) -> None:
    """Raise an error for a missing recipe.

    Args:
        recipe: Recipe name.

    Raises:
        MissingRecipeError: Always raised.
    """
    raise MissingRecipeError(
        f"Recipe '{recipe}' not found. Available recipes: "
        f"{', '.join(list(asset_build_packages.keys()))}"
    )


def _skip_lock(skip_arg: bool, cfg: str) -> bool:
    """Determine whether to skip the config file read lock.

    If config read lock skip was not forced, check if dir is writable
    and set the default to the result.

    Args:
        skip_arg: Whether skip was requested on the CLI.
        cfg: Path to the config file.

    Returns:
        Whether to skip the file lock for read.
    """
    return is_writable(os.path.dirname(cfg)) if not skip_arg else True


def make_sure_path_exists(path: str) -> None:
    """Create all directories in a path if they do not exist.

    Args:
        path: Path to create.

    Raises:
        OSError: If the path creation fails for a reason other than
            the directory already existing.
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
