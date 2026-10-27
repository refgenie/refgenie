"""CLI parsing machinery: the ``CliList`` workaround and the help rendering.

Command models and handlers live with their own family (``database.py``,
``lookup.py``, ``pull.py``, ...); what they share at handler time is in
``helpers.py``. Nothing here knows about a specific command.
"""

import argparse
import json
from typing import Annotated

from pydantic import BeforeValidator
from pydantic_settings import CliSettingsSource

from refgenie.cli.messages import COMMAND_GROUPS


# --- CliList workaround ---
#
# pydantic-settings binds a repeated CLI option to a list, but not a single
# comma-separated token. CliList adds a BeforeValidator so `-g a,b` and a JSON
# array both deserialize to a list.


def deserialize_cli_list(v):
    if isinstance(v, list):
        return v
    if isinstance(v, str):
        try:
            parsed = json.loads(v)
            if isinstance(parsed, list):
                return parsed
        except json.JSONDecodeError:
            pass
        return [x.strip() for x in v.split(",") if x.strip()]
    return v


CliList = Annotated[list, BeforeValidator(deserialize_cli_list)]


# --- CleanHelpCliSource ---
#
# Overrides help rendering: groups the root parser's subcommands (see
# COMMAND_GROUPS) and hides the auto-generated `--no-*` negation flags that
# cli_implicit_flags=True creates for every boolean field.


def grouped_format_help(parser) -> str:
    """Render grouped help output for the root parser."""
    # Find the subparsers action
    subparsers_action = None
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            subparsers_action = action
            break

    # Build a map from command name -> help text using _choices_actions
    cmd_help: dict[str, str] = {}
    if subparsers_action is not None:
        for choice_action in subparsers_action._choices_actions:
            cmd_help[choice_action.metavar] = choice_action.help or ""

    lines = []

    # Usage line
    lines.append(parser.format_usage().rstrip())

    # Description
    if parser.description:
        lines.append("")
        lines.append(parser.description)

    # Options section (just -h/--help for the root parser)
    lines.append("")
    lines.append("options:")
    for action in parser._actions:
        if not isinstance(action, argparse._SubParsersAction) and action.option_strings:
            opts = ", ".join(action.option_strings)
            help_text = action.help or ""
            lines.append(f"  {opts:<20s}{help_text}")

    # Command groups
    for group_name, commands in COMMAND_GROUPS.items():
        lines.append("")
        lines.append(f"{group_name}:")
        for cmd in commands:
            help_text = cmd_help.get(cmd, "")
            lines.append(f"    {cmd:<20s}{help_text}")

    lines.append("")
    return "\n".join(lines)


class CleanHelpCliSource(CliSettingsSource):
    def _connect_root_parser(self, *args, **kwargs):
        super()._connect_root_parser(*args, **kwargs)
        self._hide_negation_flags(self._root_parser)
        self._root_parser.format_help = lambda: grouped_format_help(self._root_parser)

    # Boolean fields whose auto-generated `--no-<leaf>` negation is a documented,
    # user-facing flag and must NOT be hidden from help. (e.g. `genome init`
    # builds the fasta asset by default; `--no-build` opts out.)
    _KEEP_NEGATION_LEAVES = frozenset({"build"})

    def _hide_negation_flags(self, parser):
        # Hide auto-generated `--no-<dest>` negation flags. pydantic-settings
        # (with cli_implicit_flags=True) pairs a `--no-<leaf>` negation with the
        # positive form for every boolean field; we drop those negations to keep
        # help output clean, EXCEPT for leaves in `_KEEP_NEGATION_LEAVES` whose
        # negation is an intentional user-facing flag (e.g. `--no-build`).
        for action in parser._actions:
            dest = getattr(action, "dest", "") or ""
            # dest may be dotted for nested fields (e.g. "genome.init.build");
            # only the leaf token contributes to the auto-negation flag name.
            leaf = dest.rsplit(".", 1)[-1].replace("_", "-")
            if leaf in self._KEEP_NEGATION_LEAVES:
                continue
            auto_negation = f"--no-{leaf}"
            kept = tuple(opt for opt in action.option_strings if opt != auto_negation)
            if kept and len(kept) != len(action.option_strings):
                action.option_strings = kept
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                for subparser in action.choices.values():
                    self._hide_negation_flags(subparser)
