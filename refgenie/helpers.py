from __future__ import annotations

import errno
import os

from refgenconf import MissingRecipeError
from ubiquerg import is_writable

from .asset_build_packages import asset_build_packages
from .exceptions import MissingFolderError


def _parse_user_kw_input(input: list[list[str]] | None) -> dict[str, str] | None:
    """Parse user input specification for build parents and input parsing.

    Args:
        input: User command line input, formatted as
            ``[[fasta=txt, test=txt], ...]``.

    Returns:
        Mapping of input names to values.
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
        MissingRecipeError: Always.
    """
    raise MissingRecipeError(
        f"Recipe '{recipe}' not found. Available recipes: "
        f"{', '.join(list(asset_build_packages.keys()))}"
    )


def _skip_lock(skip_arg: bool, cfg: str) -> bool:
    """Determine whether to skip the file lock for reading.

    If config read lock skip was not forced, checks if the directory is
    writable and uses that as the default.

    Args:
        skip_arg: Whether skip was selected on the CLI.
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
        OSError: If the path creation attempt hits an error with a code
            indicating a cause other than pre-existence.
    """
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def is_file_gzipped(path: str) -> bool:
    """Check if a file is gzipped.

    Args:
        path: Path to the file to check.

    Returns:
        Whether the file is gzipped.
    """
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"
