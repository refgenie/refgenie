"""The commands that change the local asset catalog -- `add`, `remove`,
`rename`: models and handlers.
"""

import sys
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg

from refgenie.cli.commands.framework import CliList
from refgenie.cli.commands.helpers import (
    _ASSET_REGISTRY_PATHS_DESCRIPTION,
    _validate_registry_paths,
)
from refgenie.logger import logger


class InsertModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description=_ASSET_REGISTRY_PATHS_DESCRIPTION,
    )
    path: str = Field(
        description="Relative local path to asset.",
        validation_alias=AliasChoices("p", "path"),
    )
    asset_class: str = Field(
        description="Name of the asset class of the asset.",
        validation_alias=AliasChoices("c", "asset-class"),
    )
    description: str | None = Field(
        default=None,
        description="Description of the asset.",
        validation_alias=AliasChoices("d", "description"),
    )
    seek_keys: CliList | None = Field(
        default=None,
        description="Non-path seek key values to attach to the asset. "
        "Format: name=value. Repeat for multiple keys.",
        validation_alias=AliasChoices("k", "seek-keys"),
    )


class RemoveModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description=_ASSET_REGISTRY_PATHS_DESCRIPTION,
    )
    force: bool = Field(
        default=False,
        description="Do not prompt before action, approve upfront.",
        validation_alias=AliasChoices("f", "force"),
    )
    aliases: bool = Field(
        default=False,
        description="Remove the genome alias if last asset for that genome is removed.",
        validation_alias=AliasChoices("a", "aliases"),
    )


class RenameModel(BaseModel):
    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description=_ASSET_REGISTRY_PATHS_DESCRIPTION,
    )
    new_asset_name: str = Field(
        description="New name for the asset.",
        validation_alias=AliasChoices("n", "new-asset-name"),
    )


def handle_add(cmd, refgenie) -> None:
    # Parse --seek-key name=value pairs
    custom_seek_keys = None
    if cmd.seek_keys:
        custom_seek_keys = {}
        for pair in cmd.seek_keys:
            name, _, value = pair.partition("=")
            if not name or not _:
                raise ValueError(f"Invalid seek key format: '{pair}'. Expected 'name=value'.")
            custom_seek_keys[name] = value

    parsed_asset_registry_paths = [
        refgenie.parse_asset_registry_path(p) for p in cmd.asset_registry_paths
    ]
    _validate_registry_paths(parsed_asset_registry_paths, cmd.asset_registry_paths)
    for parsed_asset_registry_path in parsed_asset_registry_paths:
        if parsed_asset_registry_path.asset is None:
            logger.error("Asset name is required. Provide it like '<genome>:<asset_group>:<asset>'")
            sys.exit(1)
        refgenie.add(
            genome_name=parsed_asset_registry_path.genome,
            asset_group_name=parsed_asset_registry_path.asset_group,
            asset_name=parsed_asset_registry_path.asset,
            path=Path(cmd.path),
            asset_class_name=cmd.asset_class,
            description=cmd.description,
            custom_seek_keys=custom_seek_keys,
        )


def handle_remove(cmd, refgenie) -> None:
    from rich.prompt import Confirm

    parsed_asset_registry_paths = [
        refgenie.parse_asset_registry_path(p) for p in cmd.asset_registry_paths
    ]
    _validate_registry_paths(parsed_asset_registry_paths, cmd.asset_registry_paths)
    touched_genomes = []
    for parsed_asset_registry_path in parsed_asset_registry_paths:
        genome_digest = refgenie.alias.resolve(parsed_asset_registry_path.genome)
        asset_name = parsed_asset_registry_path.asset or refgenie.asset.get_default(
            genome_name=parsed_asset_registry_path.genome,
            asset_group_name=parsed_asset_registry_path.asset_group,
        )
        if not cmd.force and not Confirm.ask(
            f"Are you sure you want to remove the '{genome_digest}/{parsed_asset_registry_path.asset_group}:{asset_name}' asset?"
        ):
            logger.info("Aborted by a user. Asset not removed")
            continue
        refgenie.asset.remove(
            genome_name=parsed_asset_registry_path.genome,
            genome_digest=genome_digest,
            asset_group_name=parsed_asset_registry_path.asset_group,
            asset_name=asset_name,
        )
        touched_genomes.append(genome_digest)

    if not cmd.aliases:
        return
    # Done after the whole loop, not per-iteration, so removing two assets of
    # the same genome in one invocation still triggers the cleanup.
    for genome_digest in dict.fromkeys(touched_genomes):
        if any(refgenie.asset.list_assets(genome_digests=[genome_digest])):
            continue
        for alias in refgenie.alias.get_for_genome(genome_digest=genome_digest):
            refgenie.alias.remove(alias)
            logger.info(f"Removed alias '{alias}': no assets left for {genome_digest}")


def handle_rename(cmd, refgenie) -> None:
    parsed_asset_registry_paths = [
        refgenie.parse_asset_registry_path(p) for p in cmd.asset_registry_paths
    ]
    _validate_registry_paths(parsed_asset_registry_paths, cmd.asset_registry_paths)
    for parsed_asset_registry_path in parsed_asset_registry_paths:
        asset_name = parsed_asset_registry_path.asset or refgenie.asset.get_default(
            genome_name=parsed_asset_registry_path.genome,
            asset_group_name=parsed_asset_registry_path.asset_group,
        )
        refgenie.asset.rename(
            genome_name=parsed_asset_registry_path.genome,
            asset_group_name=parsed_asset_registry_path.asset_group,
            asset_name=asset_name,
            new_asset_name=cmd.new_asset_name,
        )
