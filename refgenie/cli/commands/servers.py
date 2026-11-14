"""Talking to other refgenie servers: `subscribe`, `unsubscribe`, `catalog-export`,
`catalog-prune`."""

from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import note
from refgenie.cli.errors import fail
from refgenie.cli.framework import CliList
from refgenie.logger import logger

_GENOME_SERVER_DESCRIPTION = "One or more URLs to add to/remove from the list of subscriptions."


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
    catalog_store_name: str | None = Field(
        default=None,
        description="Name of the public store holding the genomes this node "
        "built itself (e.g. refgenie-main). Those genomes have no store; the "
        "export assigns them this one and publishes it.",
        validation_alias=AliasChoices("catalog-store-name"),
    )
    catalog_store_url: str | None = Field(
        default=None,
        description="Public root URL of the catalog store. Required with --catalog-store-name.",
        validation_alias=AliasChoices("catalog-store-url"),
    )
    catalog_store_priority: int = Field(
        default=30,
        description="Priority of the catalog store (lower wins alias ties).",
        validation_alias=AliasChoices("catalog-store-priority"),
    )


class CatalogPruneModel(BaseModel):
    source: str = Field(
        default=None,
        description="Publish-catalog artifact (path or http(s) URL) the server "
        "last imported. Defaults to $REFGENIE_CATALOG_URL.",
        validation_alias=AliasChoices("source"),
    )
    dry_run: bool = Field(
        default=False,
        description="List what would be deleted without deleting it.",
        validation_alias=AliasChoices("dry-run"),
    )


def handle_subscribe(cmd, refgenie) -> None:
    if not cmd.genome_server:
        fail("--genome-server is required: provide one or more URLs")
    refgenie.servers.subscribe(
        server_urls=cmd.genome_server,
        reset=cmd.reset,
    )


def handle_unsubscribe(cmd, refgenie) -> None:
    if not cmd.genome_server:
        fail("--genome-server is required: provide one or more URLs")
    refgenie.servers.unsubscribe(server_urls=cmd.genome_server)


def handle_catalog_export(cmd, refgenie) -> None:
    if not cmd.dest:
        fail("--dest is required: path for the publish-catalog artifact")
    if not cmd.https_prefix:
        fail("--https-prefix is required: public https base URL for the pushed assets")
    if bool(cmd.catalog_store_name) != bool(cmd.catalog_store_url):
        fail("--catalog-store-name and --catalog-store-url must be given together")
    from refgenie.publish_catalog import export_publish_catalog

    try:
        summary = export_publish_catalog(
            refgenie,
            cmd.dest,
            cmd.https_prefix,
            catalog_store_name=cmd.catalog_store_name,
            catalog_store_url=cmd.catalog_store_url,
            catalog_store_priority=cmd.catalog_store_priority,
        )
    except (ValueError, RuntimeError) as e:
        fail(f"Catalog export failed: {e}")
    note(f"Exported publish catalog to {cmd.dest}")
    for name, count in summary.items():
        logger.info(f"  {name}: {count}")


def handle_catalog_prune(cmd, refgenie) -> None:
    import os

    source = cmd.source or os.environ.get("REFGENIE_CATALOG_URL")
    if not source:
        fail("--source is required (or set $REFGENIE_CATALOG_URL): the imported catalog")
    from refgenie.publish_catalog import prune_publish_catalog

    removed = prune_publish_catalog(refgenie.database_engine, source, dry_run=cmd.dry_run)
    verb = "Would delete" if cmd.dry_run else "Deleted"
    for kind in ("asset", "assetgroup", "genome", "alias", "store"):
        note(f"{verb} {len(removed[kind])} {kind} row(s)")
        for item in removed[kind]:
            note(f"  {item}")
    if removed["kept_genome"]:
        note(f"Kept {len(removed['kept_genome'])} genome(s) that still have assets")
        for item in removed["kept_genome"]:
            note(f"  {item}")
