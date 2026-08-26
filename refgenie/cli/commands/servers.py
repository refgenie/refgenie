"""Talking to other refgenie servers: `subscribe`, `unsubscribe`, `catalog-export`."""

from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.framework import CliList
from refgenie.cli.errors import fail
from refgenie.logger import logger

_GENOME_SERVER_DESCRIPTION = (
    "One or more URLs to add to/remove from the list of subscriptions."
)


class SubscribeModel(BaseModel):
    genome_server: CliList = Field(
        default=[],
        description=_GENOME_SERVER_DESCRIPTION,
        validation_alias=AliasChoices("s", "genome-server"),
    )
    reset: bool = Field(
        default=False,
        description="Overwrite the current list of server URLs.",
        validation_alias=AliasChoices("r", "reset"),
    )


class UnsubscribeModel(BaseModel):
    genome_server: CliList = Field(
        default=[],
        description=_GENOME_SERVER_DESCRIPTION,
        validation_alias=AliasChoices("s", "genome-server"),
    )


class CatalogExportModel(BaseModel):
    dest: Path = Field(
        default=None,
        description="Path to write the publish-catalog SQLite artifact to.",
        validation_alias=AliasChoices("dest"),
    )
    https_prefix: str = Field(
        default=None,
        description="Public https base URL that mirrors the stage folder "
        "(where `refgenie push` uploaded the assets), e.g. "
        "https://<bucket>.s3.amazonaws.com/assets. Download links are "
        "served from here.",
        validation_alias=AliasChoices("https-prefix"),
    )


def handle_subscribe(cmd, refgenie) -> None:
    if not cmd.genome_server:
        fail("--genome-server is required: provide one or more URLs")
    refgenie.configuration.subscribe(
        server_urls=cmd.genome_server,
        reset=cmd.reset,
    )


def handle_unsubscribe(cmd, refgenie) -> None:
    if not cmd.genome_server:
        fail("--genome-server is required: provide one or more URLs")
    refgenie.configuration.unsubscribe(server_urls=cmd.genome_server)


def handle_catalog_export(cmd, refgenie) -> None:
    if not cmd.dest:
        fail("--dest is required: path for the publish-catalog artifact")
    if not cmd.https_prefix:
        fail("--https-prefix is required: public https base URL for the pushed assets")
    from refgenie.catalog_transfer import export_publish_catalog

    summary = export_publish_catalog(refgenie, cmd.dest, cmd.https_prefix)
    logger.info(f"Exported publish catalog to {cmd.dest}")
    for name, count in summary.items():
        logger.info(f"  {name}: {count}")
