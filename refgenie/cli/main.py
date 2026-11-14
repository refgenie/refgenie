"""Pydantic-settings CLI entry point for refgenie."""

import sys

import pydantic
from pydantic_settings import get_subcommand


def _print_validation_errors(exc: pydantic.ValidationError) -> None:
    """Format Pydantic validation errors as user-friendly CLI messages."""
    for error in exc.errors():
        loc = error.get("loc", ())
        field_name = loc[-1] if loc else "unknown"
        # Convert field names to CLI flags: single-char uses -, longer uses --
        if len(field_name) == 1:
            flag = f"-{field_name}"
        else:
            flag = f"--{field_name.replace('_', '-')}"
        error_type = error.get("type", "")
        msg = error.get("msg", "")

        if error_type == "missing":
            print(f"Error: Missing required argument '{flag}'", file=sys.stderr)
        else:
            print(f"Error: {msg} (for '{flag}')", file=sys.stderr)

    # Hint the user to run --help
    if exc.errors():
        loc = exc.errors()[0].get("loc", ())
        subcmd_parts = [p for p in loc[:-1] if isinstance(p, str)]
        if subcmd_parts:
            help_cmd = f"refgenie {' '.join(subcmd_parts)} --help"
            print(f"\nRun '{help_cmd}' for usage information.", file=sys.stderr)


def _make_cli_source(settings_cls, cli_parse_args):
    """Create a CleanHelpCliSource that hides --no-* negation flags."""
    from refgenie.cli.framework import CleanHelpCliSource

    return CleanHelpCliSource(settings_cls, cli_parse_args=cli_parse_args)


def main(test_args: list[str] | None = None) -> None:
    """Main CLI entry point.

    Args:
        test_args: If provided, parse these instead of sys.argv (for testing).
    """
    # --version / -V fast-path (before heavy imports / pydantic-settings parsing)
    if test_args is None:
        argv = sys.argv[1:]
    else:
        argv = list(test_args)
    # Only as the first argument -- `refgenie seek --version` is not a version query.
    if argv[:1] and argv[0] in ("--version", "-V"):
        from importlib.metadata import version

        try:
            print(f"refgenie {version('refgenie')}")
        except Exception:
            print("refgenie (version unknown)")
        sys.exit(0)

    # Lazy import to keep --help fast
    from refgenie.cli.commands.database import InitModel
    from refgenie.cli.dispatch import get_dispatch
    from refgenie.cli.parser import TopLevelParser

    cli_parse_args = argv if test_args is not None else True
    cli_source = _make_cli_source(TopLevelParser, cli_parse_args)

    try:
        args = TopLevelParser(_cli_settings_source=cli_source)
    except SystemExit:
        raise
    except pydantic.ValidationError as e:
        _print_validation_errors(e)
        sys.exit(1)

    subcmd = get_subcommand(args, is_required=True)

    # Lazy import heavy dependencies only when actually running commands
    from refgenie.core import Refgenie
    from refgenie.exceptions import RefgenieError
    from refgenie.logger import logger

    logger.debug(f"CLI subcmd: {type(subcmd).__name__}")
    suppress_migrations = isinstance(subcmd, InitModel)
    refgenie = Refgenie(suppress_migrations=suppress_migrations)

    dispatch = get_dispatch()
    handler = dispatch.get(type(subcmd))
    if handler is None:
        logger.error(f"Unknown command: {type(subcmd).__name__}")
        sys.exit(1)

    try:
        # One command, one post_update, however many operations it runs.
        with refgenie.batch_updates():
            handler(subcmd, refgenie)
    except (ValueError, FileNotFoundError, OSError) as e:
        logger.error(str(e))
        sys.exit(1)
    except RefgenieError as e:
        logger.error(str(e))
        sys.exit(1)
    except Exception as e:
        logger.error(f"{type(e).__name__}: {e}")
        sys.exit(1)


def main_cli() -> None:
    """Console script entry point.

    The only place that turns confirmation prompts on; elsewhere the default
    confirmer refuses rather than reading stdin (see `refgenie.utils.prompt`).
    """
    from refgenie.utils.prompt import enable_interactive_prompts

    enable_interactive_prompts()
    main()
