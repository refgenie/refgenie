"""The `list` and `listr` commands: models and handlers."""

import json

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import (
    print_table,
    _APPEND_SERVER_DESCRIPTION,
    _GENOME_DESCRIPTION,
    _GENOME_SERVER_TRANSIENT_DESCRIPTION,
    genome_digests_from,
    genome_ref,
    locate_genome,
    resolve_transient_servers,
)
from refgenie.cli.framework import CliList


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
    as_json: bool = Field(
        default=False,
        description="Print JSON instead of a table.",
        validation_alias=AliasChoices("json"),
    )


class ListrModel(BaseModel):
    genome: CliList | None = Field(
        default=None,
        description="Genome alias, e.g. hg38. A name not registered here is looked up on the servers.",
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
    as_json: bool = Field(
        default=False,
        description="Print JSON instead of a table.",
        validation_alias=AliasChoices("json"),
    )


def _listing_json(source: str, asset_data: dict, aliases_data: dict) -> dict:
    return {
        "source": source,
        "genomes": [
            {"digest": d, "aliases": aliases_data.get(d, []), "assets": assets}
            for d, assets in asset_data.items()
        ],
    }


def handle_list(cmd, refgenie) -> None:
    genome_digests = genome_digests_from(refgenie, cmd.genome, cmd.genome_digest)
    if cmd.as_json:
        listing = _listing_json("local", *refgenie.asset.list_all(genome_digests=genome_digests))
        print(json.dumps([listing], indent=2))
        return
    for table in refgenie.asset.table(genome_digests=genome_digests):
        print_table(table)


def handle_listr(cmd, refgenie) -> None:
    server_urls = resolve_transient_servers(cmd, refgenie)
    genome_digests = None
    if cmd.genome or cmd.genome_digest:
        genome_digests = [
            locate_genome(refgenie, genome_ref(a, None), server_urls)[0] for a in cmd.genome or []
        ] + [genome_ref(None, d) for d in cmd.genome_digest or []]
    if cmd.as_json:
        assets, aliases = refgenie.servers.list_assets(
            genome_digests=genome_digests, server_urls=server_urls
        )
        listings = [_listing_json(url, data, aliases.get(url, {})) for url, data in assets.items()]
        print(json.dumps(listings, indent=2))
        return
    for table in refgenie.servers.assets_table(
        genome_digests=genome_digests, server_urls=server_urls
    ):
        print_table(table)
