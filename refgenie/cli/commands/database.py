"""The `init` and `purge` commands: models and handlers."""

from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field


class InitModel(BaseModel):
    genome_folder: Path | None = Field(
        default=None,
        description="Folder for the assets refgenie manages. Unset means "
        "$REFGENIE_HOME_PATH/genomes, or the folder already configured.",
        validation_alias=AliasChoices("f", "genome-folder"),
    )
    genome_stage_folder: Path | None = Field(
        default=None,
        description="Folder where built assets are staged for serving. Unset means "
        "$REFGENIE_HOME_PATH/archives, or the folder already configured.",
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

    from refgenie.cli.commands.helpers import note

    ok = refgenie.database.init(
        genome_folder=cmd.genome_folder,
        genome_stage_folder=cmd.genome_stage_folder,
        config_version=cmd.config_version or CURRENT_CONFIG_VERSION,
        # The flag only opts out; otherwise REFGENIE_SUBSCRIBE_DEFAULT_SERVER decides.
        subscribe_default_server=None if cmd.default_server else False,
    )
    if not ok:
        return
    servers = refgenie.servers.subscriptions()
    note(f"refgenie is ready. Assets go in {refgenie.genome_folder}")
    note(
        f"Subscribed to {', '.join(servers)}. Next: `refgenie listr` shows what you can pull."
        if servers
        else "No servers subscribed. Add one with `refgenie subscribe <url>`."
    )


def handle_purge(cmd, refgenie) -> None:
    refgenie.database.purge(force=cmd.force)
