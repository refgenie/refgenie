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
    default_server: bool = Field(
        default=True,
        description="Subscribe a new configuration to the default server "
        "(https://api.refgenie.org). Use --no-default-server to start with no subscriptions.",
    )


class PurgeModel(BaseModel):
    force: bool = Field(
        default=False,
        description="Do not prompt before action, approve upfront.",
        validation_alias=AliasChoices("f", "force"),
    )


def handle_init(cmd, refgenie) -> None:
    from refgenie.const import CURRENT_CONFIG_VERSION

    refgenie.database.init(
        genome_folder=cmd.genome_folder,
        genome_stage_folder=cmd.genome_stage_folder,
        config_version=cmd.config_version or CURRENT_CONFIG_VERSION,
        # The flag only opts out; otherwise REFGENIE_SUBSCRIBE_DEFAULT_SERVER decides.
        subscribe_default_server=None if cmd.default_server else False,
    )


def handle_purge(cmd, refgenie) -> None:
    refgenie.database.purge(force=cmd.force)
