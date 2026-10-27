"""Exit codes and the single failure primitive for CLI handlers.

The error contract for every handler in ``refgenie.cli.commands``:

    Handlers raise domain exceptions or call ``fail()``. They never
    log-and-return.

A handler that logs an error and falls off the end of the function exits 0,
which tells CI, Snakemake wrappers, and shell ``&&`` chains that a failed
command succeeded. Use ``fail()`` instead.
"""

import sys
from typing import NoReturn

# Exit codes. EXIT_NOT_FOUND / EXIT_INVALID_INPUT let lookup commands
# (`seek`, `seekr`) distinguish "you asked for something that isn't there"
# from "you asked for something malformed".
EXIT_OK = 0
EXIT_GENERAL_ERROR = 1
EXIT_NOT_FOUND = 2  # missing genome / asset / seek key
EXIT_INVALID_INPUT = 3  # malformed registry path


def fail(msg: str, code: int = EXIT_GENERAL_ERROR) -> NoReturn:
    """Log ``msg`` as an error and terminate with a non-zero exit code."""
    from refgenie.logger import logger

    logger.error(msg)
    sys.exit(code)


def run_over_paths(paths, fn) -> None:
    """Apply ``fn`` to each path, accumulating the worst exit code.

    ``fn(path)`` returns an exit code (``EXIT_OK`` on success). Every path is
    attempted -- one bad path does not abort the rest -- but if any path
    failed the process exits with the highest code seen.
    """
    had_error = EXIT_OK
    for path in paths:
        had_error = max(had_error, fn(path))
    if had_error:
        sys.exit(had_error)
