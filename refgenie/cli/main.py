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
    from refgenie.cli.commands.framework import CleanHelpCliSource

    return CleanHelpCliSource(settings_cls, cli_parse_args=cli_parse_args)


def _value_consuming_flags(model) -> set[str]:
    """Derive the set of CLI flags that consume a following value from ``model``.

    Every non-``bool`` field contributes its long form (``--field-name``) plus
    every string in its ``validation_alias`` (one-character aliases become
    ``-x``, longer ones ``--xxx``). Booleans are excluded: with
    ``cli_implicit_flags=True`` they never take a value.

    Deriving this rather than hardcoding a list is the point: a value-taking
    option added to ``pull`` is picked up automatically instead of silently
    turning its value into an asset registry path.
    """
    from pydantic import AliasChoices, AliasPath

    flags: set[str] = set()
    for name, field in model.model_fields.items():
        if field.annotation is bool:
            continue
        names = {name.replace("_", "-")}
        alias = field.validation_alias
        if isinstance(alias, str):
            names.add(alias)
        elif isinstance(alias, AliasChoices):
            for choice in alias.choices:
                if isinstance(choice, str):
                    names.add(choice)
                elif isinstance(choice, AliasPath) and choice.path:
                    first = choice.path[0]
                    if isinstance(first, str):
                        names.add(first)
        if isinstance(field.alias, str):
            names.add(field.alias)
        for n in names:
            flags.add(f"-{n}" if len(n) == 1 else f"--{n}")
    return flags


def _preprocess_argv(argv: list[str]) -> list[str]:
    """Rewrite argv tokens for ergonomic invocations.

    Two transformations:
    - `data_channel` -> `data-channel` in the subcommand position only (an
      underscore alias for the kebab subcommand). Later tokens are left alone:
      a data channel may legitimately be *named* `data_channel`.
    - For `pull <paths...>`, insert `--asset-registry-paths` before the first
      positional path so pydantic-settings can bind it to the named field.
      (pydantic-settings 2.x rejects positional args with default values, so
      `asset_registry_paths` for `pull` is non-positional internally.)
    """
    if not argv:
        return argv
    out = list(argv)

    # Underscore -> kebab alias for the data_channel subcommand (argv[0] only).
    if out[0] == "data_channel":
        out[0] = "data-channel"

    # Pull positional paths -> --asset-registry-paths
    # Find the first non-flag token; if it's `pull`, collect subsequent
    # positional tokens (those not starting with `-` and not values of preceding flags)
    # and prepend `--asset-registry-paths`.
    try:
        pull_idx = out.index("pull")
    except ValueError:
        return out
    # Only treat as the subcommand if no earlier positional non-flag exists.
    for tok in out[:pull_idx]:
        if not tok.startswith("-"):
            return out
    # Walk tokens after `pull`, classify each as flag/value/positional. The set
    # of value-consuming flags is derived from PullModel, never hardcoded.
    from refgenie.cli.commands.pull import PullModel

    value_consuming = _value_consuming_flags(PullModel)
    new_tail: list[str] = []
    paths: list[str] = []
    i = pull_idx + 1
    skip_next = False
    while i < len(out):
        tok = out[i]
        if skip_next:
            new_tail.append(tok)
            skip_next = False
        elif tok.startswith("-"):
            new_tail.append(tok)
            # Detect value-consuming long form like --genome=value (no extra arg)
            base = tok.split("=", 1)[0]
            if "=" not in tok and base in value_consuming:
                skip_next = True
        else:
            paths.append(tok)
        i += 1
    if paths and "--asset-registry-paths" not in new_tail:
        new_tail = ["--asset-registry-paths", ",".join(paths)] + new_tail
    return out[: pull_idx + 1] + new_tail


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

    # Argv preprocessing: pull positional shim + data_channel underscore alias
    argv = _preprocess_argv(argv)

    # Lazy import to keep --help fast
    from refgenie.cli.commands.database import InitModel
    from refgenie.cli.dispatch import get_dispatch
    from refgenie.cli.parser import TopLevelParser

    cli_parse_args = argv if (test_args is not None or argv != sys.argv[1:]) else True
    cli_source = _make_cli_source(TopLevelParser, cli_parse_args)

    try:
        args = TopLevelParser(_cli_settings_source=cli_source)
    except SystemExit:
        raise  # Let argparse errors propagate normally
    except pydantic.ValidationError as e:
        _print_validation_errors(e)
        sys.exit(1)

    subcmd = get_subcommand(args, is_required=True)

    # Lazy import heavy dependencies only when actually running commands
    from refgenie.exceptions import RefgenieError
    from refgenie.logger import logger
    from refgenie.core import Refgenie

    logger.debug(f"CLI subcmd: {type(subcmd).__name__}")
    suppress_migrations = isinstance(subcmd, InitModel)
    refgenie = Refgenie(suppress_migrations=suppress_migrations)

    dispatch = get_dispatch()
    handler = dispatch.get(type(subcmd))
    if handler is None:
        logger.error(f"Unknown command: {type(subcmd).__name__}")
        sys.exit(1)

    try:
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

    This is the one context in refgenie that owns a terminal, so it is the one
    place allowed to turn confirmation prompts on. Everywhere else -- library
    callers, the server, the MCP tools, tests -- the default confirmer refuses
    rather than reading stdin. See `refgenie.utils.prompt`.
    """
    from refgenie.utils.prompt import enable_interactive_prompts

    enable_interactive_prompts()
    main()
