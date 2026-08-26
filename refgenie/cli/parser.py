"""The top-level CLI parser, assembling every command family's models."""

from pydantic import Field
from pydantic_settings import BaseSettings, CliSubCommand, SettingsConfigDict

from refgenie.cli.commands.alias import AliasParser
from refgenie.cli.commands.asset import AssetParser
from refgenie.cli.commands.asset_class import AssetClassParser
from refgenie.cli.commands.build import BuildModel
from refgenie.cli.commands.config import ConfigParser
from refgenie.cli.commands.data_channel import DataChannelParser
from refgenie.cli.commands.curate import InsertModel, RemoveModel, RenameModel
from refgenie.cli.commands.database import InitModel, PurgeModel
from refgenie.cli.commands.generate import GenerateParser
from refgenie.cli.commands.genome import GenomeParser
from refgenie.cli.commands.getseq import GetseqModel
from refgenie.cli.commands.listing import ListModel, ListrModel
from refgenie.cli.commands.lookup import CompareModel, IdModel, SeekModel, SeekrModel
from refgenie.cli.commands.populate import PopulateModel, PopulaterModel
from refgenie.cli.commands.pull import MirrorModel, PullModel
from refgenie.cli.commands.push import PushModel
from refgenie.cli.commands.recipe import RecipeParser
from refgenie.cli.commands.remote import RemoteParser
from refgenie.cli.commands.store import StoreParser
from refgenie.cli.commands.serve import DashModel, ServeModel
from refgenie.cli.commands.stage import StageParser
from refgenie.cli.commands.servers import (
    CatalogExportModel,
    SubscribeModel,
    UnsubscribeModel,
)
from refgenie.cli.messages import SUBPARSER_MESSAGES


class TopLevelParser(BaseSettings):
    """refgenie - reference genome asset manager"""

    model_config = SettingsConfigDict(
        cli_parse_args=True,
        cli_prog_name="refgenie",
        cli_kebab_case=True,
        cli_implicit_flags=True,
        cli_hide_none_type=True,
    )

    # Config/subscription and catalog commands
    init: CliSubCommand[InitModel] = Field(description=SUBPARSER_MESSAGES["init"])
    purge: CliSubCommand[PurgeModel] = Field(description=SUBPARSER_MESSAGES["purge"])
    list: CliSubCommand[ListModel] = Field(description=SUBPARSER_MESSAGES["list"])
    subscribe: CliSubCommand[SubscribeModel] = Field(description=SUBPARSER_MESSAGES["subscribe"])
    unsubscribe: CliSubCommand[UnsubscribeModel] = Field(
        description=SUBPARSER_MESSAGES["unsubscribe"]
    )
    catalog_export: CliSubCommand[CatalogExportModel] = Field(
        alias="catalog-export",
        description=SUBPARSER_MESSAGES["catalog_export"],
    )

    # Asset lookup and transfer commands
    seek: CliSubCommand[SeekModel] = Field(description=SUBPARSER_MESSAGES["seek"])
    seekr: CliSubCommand[SeekrModel] = Field(description=SUBPARSER_MESSAGES["seekr"])
    remove: CliSubCommand[RemoveModel] = Field(description=SUBPARSER_MESSAGES["remove"])
    rename: CliSubCommand[RenameModel] = Field(description=SUBPARSER_MESSAGES["rename"])
    id: CliSubCommand[IdModel] = Field(description=SUBPARSER_MESSAGES["id"])
    add: CliSubCommand[InsertModel] = Field(description=SUBPARSER_MESSAGES["add"])
    getseq: CliSubCommand[GetseqModel] = Field(description=SUBPARSER_MESSAGES["getseq"])
    pull: CliSubCommand[PullModel] = Field(description=SUBPARSER_MESSAGES["pull"])
    listr: CliSubCommand[ListrModel] = Field(description=SUBPARSER_MESSAGES["listr"])
    compare: CliSubCommand[CompareModel] = Field(description=SUBPARSER_MESSAGES["compare"])
    populate: CliSubCommand[PopulateModel] = Field(description=SUBPARSER_MESSAGES["populate"])
    populater: CliSubCommand[PopulaterModel] = Field(description=SUBPARSER_MESSAGES["populater"])
    mirror: CliSubCommand[MirrorModel] = Field(description=SUBPARSER_MESSAGES["mirror"])

    # Build command
    build: CliSubCommand[BuildModel] = Field(description=SUBPARSER_MESSAGES["build"])

    # Serve and Dash commands
    serve: CliSubCommand[ServeModel] = Field(description=SUBPARSER_MESSAGES["serve"])
    dash: CliSubCommand[DashModel] = Field(description=SUBPARSER_MESSAGES["dash"])

    # Push command
    push: CliSubCommand[PushModel] = Field(description=SUBPARSER_MESSAGES["push"])

    # Nested command groups
    alias: CliSubCommand[AliasParser] = Field(description=SUBPARSER_MESSAGES["alias"])
    config: CliSubCommand[ConfigParser] = Field(description=SUBPARSER_MESSAGES["config"])
    recipe: CliSubCommand[RecipeParser] = Field(description=SUBPARSER_MESSAGES["recipe"])
    asset_class: CliSubCommand[AssetClassParser] = Field(
        alias="asset-class",
        description=SUBPARSER_MESSAGES["asset_class"],
    )
    stage: CliSubCommand[StageParser] = Field(description=SUBPARSER_MESSAGES["stage"])
    data_channel: CliSubCommand[DataChannelParser] = Field(
        alias="data-channel",
        description=SUBPARSER_MESSAGES["data_channel"],
    )
    generate: CliSubCommand[GenerateParser] = Field(description=SUBPARSER_MESSAGES["generate"])
    remote: CliSubCommand[RemoteParser] = Field(description=SUBPARSER_MESSAGES["remote"])
    store: CliSubCommand[StoreParser] = Field(description=SUBPARSER_MESSAGES["store"])
    genome: CliSubCommand[GenomeParser] = Field(description=SUBPARSER_MESSAGES["genome"])
    asset: CliSubCommand[AssetParser] = Field(
        description="Asset operations (alias group: provides 'asset list')."
    )
