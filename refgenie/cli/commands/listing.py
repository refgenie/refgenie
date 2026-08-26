"""The `list` and `listr` commands: models and handlers."""

from pydantic import AliasChoices, BaseModel, Field
from rich import print as rprint

from refgenie.cli.commands.framework import CliList
from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _GENOME_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    _resolve_genome_arg,
    resolve_transient_servers,
)


class ListModel(BaseModel):
    genome: CliList | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION,
        validation_alias=AliasChoices("g", "genome"),
    )


class ListrModel(BaseModel):
    genome: CliList | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION,
        validation_alias=AliasChoices("g", "genome"),
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


def handle_list(cmd, refgenie) -> None:
    for table in refgenie.asset.table(genome_names=cmd.genome):
        rprint(table)


def handle_listr(cmd, refgenie) -> None:
    genome_digests = None
    if cmd.genome:
        genome_digests = [_resolve_genome_arg(refgenie, g) for g in cmd.genome]
    for table in refgenie.asset.remote_table(
        genome_digests=genome_digests,
        server_urls=resolve_transient_servers(cmd, refgenie),
    ):
        rprint(table)
