"""The `build` command: model and handler."""

import sys

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliPositionalArg
from rich import print as rprint

from refgenie.cli.commands.framework import CliList
from refgenie.cli.commands.helpers import (
    _ASSET_REGISTRY_PATHS_DESCRIPTION,
    _validate_registry_paths,
    data_channel_hint,
    no_data_channel_synced,
)
from refgenie.logger import logger


class BuildModel(BaseModel):
    """Model for the 'build' subcommand arguments."""

    asset_registry_paths: CliPositionalArg[list[str]] = Field(
        description=_ASSET_REGISTRY_PATHS_DESCRIPTION,
    )

    asset_description: str | None = Field(
        default=None,
        description="Add asset level description (e.g. built with version 0.3.2).",
        validation_alias=AliasChoices("asset-description"),
    )

    recipe_name: str | None = Field(
        default=None,
        description="Provide a recipe to use.",
        validation_alias=AliasChoices("recipe-name"),
    )

    recipe_version: str | None = Field(
        default=None,
        description="Provide a recipe version to use.",
        validation_alias=AliasChoices("recipe-version"),
    )

    docker: bool = Field(
        default=False,
        description="Run all commands in the refgenie docker container.",
        validation_alias=AliasChoices("d", "docker"),
    )

    pull_parents: bool = Field(
        default=False,
        description="Automatically pull the default parent asset if required but not provided.",
        validation_alias=AliasChoices("pull-parents"),
    )

    requirements: bool = Field(
        default=False,
        description="Show the build requirements for the specified asset and exit.",
        validation_alias=AliasChoices("q", "requirements"),
    )

    stage: bool = Field(
        default=False,
        description="Stage the asset after building. Requires the genome stage folder to be set.",
    )

    push_to: list[str] | None = Field(
        default=None,
        description="Remote names/IDs to create push intent records for after staging.",
        validation_alias=AliasChoices("push-to"),
    )

    pipeline_kwargs: CliList | None = Field(
        default=None,
        description="Extra arguments to pass to the build pipeline. Format: arg_name=arg_val arg_name1=arg_val1",
        validation_alias=AliasChoices("pipeline-kwargs"),
    )

    assets: CliList | None = Field(
        default=None,
        description="Override the default genome, asset and tag of the parents "
        "(e.g. fasta=hg38/fasta:default gtf=mm10/gencode_gtf:default).",
    )

    files: CliList | None = Field(
        default=None,
        description="Provide paths to the required files (e.g. fasta=/path/to/file.fa.gz).",
    )

    params: CliList | None = Field(
        default=None,
        description="Provide required parameter values (e.g. param1=value1).",
    )

    volumes: CliList | None = Field(
        default=None,
        description="If using docker, also mount these folders as volumes.",
    )


def handle_build(cmd, refgenie) -> None:
    from refgenie.exceptions import MissingAliasError, MissingAssetClassError, MissingRecipeError
    from refgenie.models import BuildParams
    from refgenie.utils.io import parse_user_kw_input

    parsed_asset_registry_paths = [
        refgenie.parse_asset_registry_path(p) for p in cmd.asset_registry_paths
    ]
    _validate_registry_paths(parsed_asset_registry_paths, cmd.asset_registry_paths)
    recipe_name = None
    if cmd.recipe_name:
        if len(parsed_asset_registry_paths) > 1:
            logger.error("Recipes cannot be specified for multi-asset builds")
            sys.exit(1)
        recipe_name = cmd.recipe_name
    if cmd.requirements:
        recipe_names = []
        for parsed_asset_registry_path in parsed_asset_registry_paths:
            recipe = refgenie.recipe.get(
                recipe_name=recipe_name or parsed_asset_registry_path.asset_group,
                recipe_version=cmd.recipe_version,
            )
            recipe_names.append(recipe.name)
        rprint(refgenie.recipe.table(recipe_names=recipe_names))
        sys.exit(0)

    pipeline_kwargs = parse_user_kw_input(cmd.pipeline_kwargs)
    specified_files = parse_user_kw_input(cmd.files)
    specified_params = parse_user_kw_input(cmd.params)
    specified_assets = parse_user_kw_input(cmd.assets)
    for parsed_asset_registry_path in parsed_asset_registry_paths:
        build_params = BuildParams(
            assets=specified_assets,
            params=specified_params,
            files=specified_files,
        )
        try:
            refgenie.build_asset(
                genome_name=parsed_asset_registry_path.genome,
                asset_group_name=parsed_asset_registry_path.asset_group,
                asset_name=parsed_asset_registry_path.asset,
                recipe_name=recipe_name or parsed_asset_registry_path.asset_group,
                recipe_version=cmd.recipe_version,
                params=build_params,
                stage=cmd.stage,
                push_to=cmd.push_to,
                asset_description=cmd.asset_description,
                pull_parents=cmd.pull_parents,
                pipeline_kwargs=pipeline_kwargs,
                docker=cmd.docker,
                docker_volumes=cmd.volumes or [],
            )
        except MissingAliasError:
            genome = parsed_asset_registry_path.genome
            if parsed_asset_registry_path.asset_group == "fasta":
                fasta_path = specified_files.get("fasta", "<file>")
                logger.info(f"Hint: refgenie genome init --fasta {fasta_path} --name {genome}")
            sys.exit(1)
        except (MissingRecipeError, MissingAssetClassError) as e:
            logger.error(str(e))
            if no_data_channel_synced(refgenie):
                genome = parsed_asset_registry_path.genome
                asset_group = parsed_asset_registry_path.asset_group
                logger.info(data_channel_hint(f"refgenie build {genome}/{asset_group}"))
            sys.exit(1)
