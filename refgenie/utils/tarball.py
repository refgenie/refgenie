import os
import shutil
from pathlib import Path
from subprocess import run

from refgenie.utils.console import CONSOLE
from refgenie.logger import logger


def tar(path: Path, output: Path, exclude_files: list[str] | None = None) -> None:
    """
    Tar the asset directory.

    Args:
        path: Path to the asset directory to tar.
        output: Path to the output tarball.
        exclude_files: Optional list of filenames to exclude (e.g. colocation symlinks).
    """
    pth, dir_name = os.path.split(path)
    if not path.exists():
        raise FileNotFoundError(f"Entity '{path}' does not exist")
    # tar gzip the asset; the directory holds asset content only, so the only
    # exclusions are the caller's (e.g. colocation symlinks).
    excludes = " ".join(f"--exclude '{f}'" for f in exclude_files or [])
    compressor = "pigz" if shutil.which("pigz") else "gzip"
    cmd = f"tar {excludes} -C {pth} -cvf - {dir_name} | {compressor} > {output}"
    logger.info(f"Running: {cmd}")
    with CONSOLE.status(f"[green]Creating tarball in {output}..."):
        # check=True: a silently truncated tarball would be digested, recorded
        # as a StagedAsset and pushed as if it were complete.
        run(cmd, shell=True, check=True)


def read_build_commands(build_dir: Path) -> list[str]:
    """
    Read build commands from the build bookkeeping directory.

    Args:
        build_dir: Path to the build directory (see ``utils.build.get_build_dir``).

    Returns:
        List of build commands. Empty if the asset was pulled rather than built,
        in which case no build directory exists.
    """
    if not build_dir.is_dir():
        logger.debug(f"No build directory found at {build_dir}. No build commands to read.")
        return []
    cmds_paths = list(build_dir.glob("*_commands.sh"))
    if not len(cmds_paths) == 1:
        logger.warning(f"Expected one build commands, got {len(cmds_paths)}")
        return []
    with open(cmds_paths[0], "r") as f:
        build_commands_list = f.read().split("\n")
        return [cmd for cmd in build_commands_list if cmd and not cmd.startswith("#")]


def copy_asset_dir(input_dir: Path, target_dir: Path) -> None:
    """
    Copy the asset directory.

    Args:
        input_dir: Path to the directory to copy the asset dir from.
        target_dir: Path to the directory to copy the asset dir to.
    """
    if not input_dir.exists():
        raise FileNotFoundError(f"Asset directory not found: {input_dir}")
    cmd = f"rsync -rvL {input_dir}/ {target_dir}/"
    logger.info(f"Running: {cmd}")
    with CONSOLE.status(f"[green]Copying asset directory from {input_dir} to {target_dir}..."):
        run(cmd, shell=True, check=True)


def get_external_symlinks(asset_dir: Path) -> list[str]:
    """Find symlinks in asset_dir that point outside the asset directory.

    These are colocation symlinks created during build that should be excluded
    from archives (they'll be recreated after pull).

    Args:
        asset_dir: Path to the asset directory.

    Returns:
        List of filenames (relative to asset_dir) that are external symlinks.
    """
    external = []
    for root, _dirs, files in os.walk(asset_dir):
        for f in files:
            file_path = Path(root) / f
            if not file_path.is_symlink():
                continue
            # Resolve the symlink target relative to its location
            target = (file_path.parent / os.readlink(file_path)).resolve()
            try:
                target.relative_to(asset_dir.resolve())
            except ValueError:
                # Target is outside asset_dir — this is a colocation symlink
                external.append(os.path.relpath(file_path, asset_dir))
    return external


def read_asset_dir_contents(target_path: Path, exclude_files: set[str] | None = None) -> list[str]:
    """
    Create a file tree with contents of the asset directory.

    The asset directory holds asset content only — build bookkeeping lives in
    a separate ``builds/`` tree — so everything found here is listed.

    Args:
        target_path: Path to the asset directory (if a directory) or a file
            whose parent directory will be scanned.
        exclude_files: Optional set of filenames to exclude from the listing
            (e.g., the tarball name).

    Returns:
        Sorted list of relative file paths in the asset directory.
    """
    if target_path.is_dir():
        asset_dir = target_path
    else:
        asset_dir = target_path.parent

    exclude = exclude_files or set()
    file_set = set()
    for dp, _dirs, fn in os.walk(asset_dir):
        for f in fn:
            file_set.add(os.path.relpath(os.path.join(dp, f), asset_dir))
    file_set -= exclude
    return sorted(list(file_set))
