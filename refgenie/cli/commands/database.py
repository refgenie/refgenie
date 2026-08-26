"""The `init` and `purge` commands: models and handlers."""

from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field


class InitModel(BaseModel):
    genome_folder: Path | None = Field(
        default=None,
        description="Absolute path to parent folder refgenie-managed assets.",
        validation_alias=AliasChoices("f", "genome-folder"),
    )
    genome_stage_folder: Path | None = Field(
        default=None,
        description="Absolute path to parent stage folder refgenie-managed assets; "
        "used by refgenieserver.",
        validation_alias=AliasChoices("a", "genome-stage-folder"),
    )
    config_version: str | None = Field(
        default=None,
        description="Config version to initialize the config file with.",
        validation_alias=AliasChoices("v", "config-version"),
    )


class PurgeModel(BaseModel):
    force: bool = Field(
        default=False,
        description="Do not prompt before action, approve upfront.",
        validation_alias=AliasChoices("f", "force"),
    )


def handle_init(cmd, refgenie) -> None:
    from refgenie.const import CURRENT_CONFIG_VERSION

    refgenie.init(
        genome_folder=cmd.genome_folder,
        genome_stage_folder=cmd.genome_stage_folder,
        config_version=cmd.config_version or CURRENT_CONFIG_VERSION,
    )


def handle_purge(cmd, refgenie) -> None:
    refgenie.purge(force=cmd.force)
