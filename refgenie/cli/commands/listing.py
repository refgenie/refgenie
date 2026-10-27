"""The `list` and `listr` commands: models and handlers."""

from pydantic import AliasChoices, BaseModel, Field
from rich import print as rprint

from refgenie.cli.framework import CliList
from refgenie.cli.commands.helpers import (
    _APPEND_SERVER_DESCRIPTION,
    _GENOME_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    genome_digests_from,
    resolve_transient_servers,
)
from refgenie.cli.errors import EXIT_NOT_FOUND, fail


class ListModel(BaseModel):
    genome: CliList | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION,
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: CliList | None = Field(
        default=None,
        description="Genome digests, used in place of aliases.",
        validation_alias=AliasChoices("genome-digest"),
    )


class ListrModel(BaseModel):
    genome: CliList | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION,
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: CliList | None = Field(
        default=None,
        description="Genome digests, used in place of aliases. The genomes need not exist locally.",
        validation_alias=AliasChoices("genome-digest"),
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
    genome_digests = genome_digests_from(refgenie, cmd.genome, cmd.genome_digest)
    for table in refgenie.asset.table(genome_digests=genome_digests):
        rprint(table)


def handle_listr(cmd, refgenie) -> None:
    from refgenie.exceptions import MissingAliasError

    try:
        genome_digests = genome_digests_from(refgenie, cmd.genome, cmd.genome_digest)
    except MissingAliasError as e:
        fail(
            f"{e} -g takes local aliases. To list a genome that has no local alias, "
            f"use --genome-digest.",
            EXIT_NOT_FOUND,
        )
    for table in refgenie.servers.assets_table(
        genome_digests=genome_digests,
        server_urls=resolve_transient_servers(cmd, refgenie),
    ):
        rprint(table)
