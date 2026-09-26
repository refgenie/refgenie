"""The `stage` command group: models and handlers."""

from collections.abc import Callable

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg, CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.commands.helpers import (
    _GENOME_DIGEST_DESCRIPTION,
    _validate_registry_paths,
    genome_digest_from,
)
from refgenie.cli.errors import fail


class StageListModel(BaseModel):
    """stage list: list staged assets."""

    pass


class StageRemoveModel(BaseModel):
    """stage remove: unstage an asset."""

    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description="Registry path strings identifying assets (e.g. hg38/fasta).",
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION,
        validation_alias=AliasChoices("genome-digest"),
    )


class StageAddModel(BaseModel):
    """stage add: stage an asset."""

    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description="Registry path strings identifying assets (e.g. hg38/fasta).",
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION,
        validation_alias=AliasChoices("genome-digest"),
    )


class StageParser(BaseModel):
    """Intermediate parser for stage add/remove/list subcommands."""

    add: CliSubCommand[StageAddModel] = Field(description="Add an asset to staging.")
    remove: CliSubCommand[StageRemoveModel] = Field(description="Remove an asset from staging.")
    list: CliSubCommand[StageListModel] = Field(description="List staged assets.")


def handle_stage_list(cmd, refgenie) -> None:
    rprint(refgenie.stage.table())


def handle_stage_add(cmd, refgenie) -> None:
    parsed_asset_registry_paths = [
        refgenie.parse_asset_registry_path(p) for p in cmd.asset_registry_paths
    ]
    _validate_registry_paths(
        parsed_asset_registry_paths, cmd.asset_registry_paths, cmd.genome_digest
    )
    for parsed_asset_registry_path in parsed_asset_registry_paths:
        genome_digest = genome_digest_from(
            refgenie, parsed_asset_registry_path.genome, cmd.genome_digest
        )
        asset_name = parsed_asset_registry_path.asset or refgenie.asset.group.get_default(
            genome_digest=genome_digest,
            asset_group_name=parsed_asset_registry_path.asset_group,
        )
        asset = refgenie.asset.get(
            genome_digest=genome_digest,
            asset_group_name=parsed_asset_registry_path.asset_group,
            asset_name=asset_name,
        )
        refgenie.stage.create(
            asset=asset,
            genome_folder=refgenie.genome_folder,
            genome_stage_folder=refgenie.genome_stage_folder,
            # The builds/ tree is keyed by the alias used at BUILD time, which
            # need not be the alias named here, so search all aliases for the
            # genome. Returns None for a pulled asset, which stages with empty
            # build_commands.
            build_dir=refgenie.asset.tree.find_build_dir(
                genome_digest=genome_digest,
                asset_group_name=parsed_asset_registry_path.asset_group,
                asset_name=asset_name,
            ),
        )


def handle_stage_remove(cmd, refgenie) -> None:
    parsed_asset_registry_paths = [
        refgenie.parse_asset_registry_path(p) for p in cmd.asset_registry_paths
    ]
    _validate_registry_paths(
        parsed_asset_registry_paths, cmd.asset_registry_paths, cmd.genome_digest
    )
    for parsed_asset_registry_path in parsed_asset_registry_paths:
        genome_digest = genome_digest_from(
            refgenie, parsed_asset_registry_path.genome, cmd.genome_digest
        )
        asset = refgenie.asset.get(
            genome_digest=genome_digest,
            asset_group_name=parsed_asset_registry_path.asset_group,
            asset_name=parsed_asset_registry_path.asset
            or refgenie.asset.group.get_default(
                genome_digest=genome_digest,
                asset_group_name=parsed_asset_registry_path.asset_group,
            ),
        )
        refgenie.stage.remove(asset_digest=asset.digest)


STAGE_DISPATCH: dict[type, Callable] = {
    StageListModel: handle_stage_list,
    StageAddModel: handle_stage_add,
    StageRemoveModel: handle_stage_remove,
}


def handle_stage_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = STAGE_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown stage subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
