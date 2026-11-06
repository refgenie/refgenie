"""Exit codes and ``fail()`` for CLI handlers.

Handlers raise domain exceptions or call ``fail()``; never log-and-return,
which exits 0 on failure.
"""

import sys
from typing import NoReturn

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

    ``fn(path)`` returns an exit code. Every path is attempted; if any failed,
    the process exits with the highest code seen.
    """
    had_error = EXIT_OK
    for path in paths:
        had_error = max(had_error, fn(path))
    if had_error:
        sys.exit(had_error)
