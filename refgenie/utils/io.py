"""
Small CLI input/output helpers: read YAML, feed files or stdin line by line to a
callback, pretty-print YAML, and parse ``key=value`` CLI arguments into typed values.
"""

import sys
from collections.abc import Callable
from pathlib import Path

from rich import print
from rich.syntax import Syntax
from ubiquerg import convert_value
from yaml import safe_load

from refgenie.logger import logger
from refgenie.utils.http import make_client


def read_yaml(path: Path | str) -> dict:
    """
    Read a YAML either from a local file or from a remote URL.

    Local paths may be given as either a ``Path`` or a plain ``str``. URLs are
    fetched with httpx, following redirects and carrying an explicit
    ``User-Agent`` header (required to get through Cloudflare, which otherwise
    403s the default UA).

    Args:
        path: The path to the YAML file (Path or str) or an http(s) URL.

    Returns:
        dict: The YAML content.
    """
    if isinstance(path, str) and path.startswith(("http://", "https://")):
        logger.info(f"Reading YAML from URL: {path}")
        with make_client() as client:
            response = client.get(path)
            response.raise_for_status()
            return safe_load(response.text)
    logger.info(f"Reading YAML from file: {path}")
    return safe_load(Path(path).read_text())


def populate_file(file_path: Path, pop_fun: Callable[[str], str]) -> None:
    """
    Read a file line by line, apply pop_fun to each line, and write the
    populated result to stdout. The input file is not modified.

    Args:
        file_path: Path to the input file containing refgenie registry paths.
        pop_fun: A function that populates refgenie registry paths in objects.
    """
    logger.debug(f"Populating file: {file_path}")
    with open(file_path) as fp:
        for line in fp:
            sys.stdout.write(pop_fun(line))


def populate_stdin(pop_fun: Callable[[str], str]) -> None:
    """
    Read lines from stdin, apply pop_fun to each, and write the result to
    stdout; a line of 'q', 'quit', or 'exit' stops processing.

    Args:
        pop_fun: A function that populates refgenie registry paths in objects.
    """
    for line in sys.stdin:
        if line.rstrip() in ["q", "quit", "exit"]:
            break
        sys.stdout.write(pop_fun(line))


def cli_show_yaml(string: str) -> None:
    """
    Print a yaml string with syntax highlighting.

    Args:
        string: The yaml string to print.
    """
    print(
        Syntax(
            string,
            "yaml",
            padding=1,
            theme="ansi_light",
            background_color="default",
        )
    )


def parse_user_kw_input(input: list[str] | None) -> dict[str, str]:
    """Parse user key=value input from CLI.

    Args:
        input: List of 'key=value' strings from CLI, or None.

    Returns:
        Dict mapping keys to values.
    """
    if not input:
        return {}
    return {k: v for k, v in (x.split("=", 1) for x in input if "=" in x)}


def coerce_cli_kwargs(kwargs: dict[str, str]) -> dict[str, str | bool | int | float]:
    """Coerce all values in a CLI kwargs dict to appropriate Python types."""
    return {k: convert_value(v) for k, v in kwargs.items()}
