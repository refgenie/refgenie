"""The `populate` and `populater` commands: models and handlers."""

from functools import partial
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    resolve_transient_servers,
)
from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.cli.framework import CliList


class PopulateModel(BaseModel):
    file: str | None = Field(
        default=None,
        description="File with registry paths to populate.",
        validation_alias=AliasChoices("f", "file"),
    )


class PopulaterModel(BaseModel):
    file: str | None = Field(
        default=None,
        description="File with registry paths to populate.",
        validation_alias=AliasChoices("f", "file"),
    )
    genome_server: CliList | None = Field(
        default=None,
        description=_GENOME_SERVER_TRANSIENT_DESCRIPTION,
        validation_alias=AliasChoices("s", "genome-server"),
    )
    append_server: bool = Field(
        default=False,
        description=_APPEND_SERVER_DESCRIPTION,
        validation_alias=AliasChoices("p", "append-server"),
    )


def _populate(file: str | None, pop_fun) -> None:
    """Populate the file (or stdin) to stdout; fail if any path was left unresolved."""
    from refgenie.utils.io import populate_file, populate_stdin

    unresolved = populate_file(Path(file), pop_fun) if file else populate_stdin(pop_fun)
    if unresolved:
        fail(f"{unresolved} line(s) kept a refgenie:// path that was not populated", EXIT_NOT_FOUND)


def handle_populate(cmd, refgenie) -> None:
    _populate(cmd.file, refgenie.populate)


def handle_populater(cmd, refgenie) -> None:
    servers = resolve_transient_servers(cmd, refgenie)
    _populate(cmd.file, partial(refgenie.populater, server_urls=servers))
