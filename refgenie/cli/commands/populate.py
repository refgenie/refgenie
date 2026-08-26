"""The `populate` and `populater` commands: models and handlers."""

from functools import partial
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.framework import CliList
from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    resolve_transient_servers,
)


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


def handle_populate(cmd, refgenie) -> None:
    from refgenie.utils.io import populate_file, populate_stdin

    file_path = Path(cmd.file) if cmd.file else None
    if file_path is not None:
        populate_file(file_path=file_path, pop_fun=refgenie.populate)
    else:
        populate_stdin(pop_fun=refgenie.populate)


def handle_populater(cmd, refgenie) -> None:
    from refgenie.utils.io import populate_file, populate_stdin

    pop_fun = partial(
        refgenie.populater,
        server_urls=resolve_transient_servers(cmd, refgenie),
    )
    file_path = Path(cmd.file) if cmd.file else None
    if file_path is not None:
        populate_file(file_path=file_path, pop_fun=pop_fun)
    else:
        populate_stdin(pop_fun=pop_fun)
