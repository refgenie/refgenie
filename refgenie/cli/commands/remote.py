"""The `remote` command group: models and handlers."""

from collections.abc import Callable
from typing import Literal

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.logger import logger


class RemoteAddModel(BaseModel):
    """remote add: add a remote."""

    type: Literal["s3", "http", "https"] = Field(description="Type of the remote.")
    prefix: str = Field(description="Prefix/identifier for the remote.")
    description: str = Field(description="Description of the remote.")
    push_command: str | None = Field(
        default=None,
        description=(
            "Shell command template for pushing assets. "
            "Placeholders: {local_path}, {relative_path}, {prefix}, {genome_stage_folder}. "
            "Example: 'aws s3 cp {local_path} s3://bucket/{relative_path}'"
        ),
        validation_alias=AliasChoices("push-command"),
    )


class RemoteRemoveModel(BaseModel):
    """remote remove: remove a remote."""

    type: Literal["s3", "http", "https"] = Field(description="Type of the remote to remove.")


class RemoteListModel(BaseModel):
    """remote list: list all remotes."""

    pass


class RemoteStatusModel(BaseModel):
    """remote status: show push status of remote-asset links."""

    remote: str | None = Field(
        default=None,
        description="Show status for only this remote (by name or id).",
        validation_alias=AliasChoices("r", "remote"),
    )


class RemoteParser(BaseModel):
    """Intermediate parser for remote subcommands."""

    add: CliSubCommand[RemoteAddModel] = Field(description="Add a remote.")
    remove: CliSubCommand[RemoteRemoveModel] = Field(description="Remove a remote.")
    list: CliSubCommand[RemoteListModel] = Field(description="List all remotes.")
    status: CliSubCommand[RemoteStatusModel] = Field(description="Show push status.")


def handle_remote_list(cmd, refgenie) -> None:
    if table := refgenie.configuration.remote_table():
        rprint(table)
    else:
        logger.info("No remotes configured.")


def handle_remote_add(cmd, refgenie) -> None:
    from refgenie.db.tables import RemoteType
    from refgenie.exceptions import RefgenieError

    try:
        refgenie.configuration.add_remote(
            remote_type=RemoteType(cmd.type),
            prefix=cmd.prefix,
            description=cmd.description,
            push_command=cmd.push_command,
        )
    except (ValueError, RefgenieError) as e:
        fail(f"Failed to add remote: {e}")


def handle_remote_remove(cmd, refgenie) -> None:
    from refgenie.db.tables import RemoteType
    from refgenie.exceptions import MissingRemoteError, RefgenieError

    try:
        refgenie.configuration.remove_remote(remote_type=RemoteType(cmd.type))
    except MissingRemoteError as e:
        fail(f"Failed to remove remote: {e}", EXIT_NOT_FOUND)
    except (ValueError, RefgenieError) as e:
        fail(f"Failed to remove remote: {e}")


def handle_remote_status(cmd, refgenie) -> None:
    """Show push status of remote-asset links."""
    status = refgenie.configuration.remote_status(remote=cmd.remote)
    if not status:
        logger.info("No remotes configured.")
        return

    for entry in status.values():
        remote, pushed, unpushed = entry["remote"], entry["pushed"], entry["unpushed"]
        total = len(pushed) + len(unpushed)

        print(f"\nRemote: {remote.description or remote.type.value} (id={remote.id})")
        print(f"  Type:     {remote.type.value}")
        print(f"  Prefix:   {remote.prefix}")
        print(f"  Total:    {total} asset links")
        if total > 0:
            print(f"  Pushed:   {len(pushed)}")
            print(f"  Unpushed: {len(unpushed)}")
        if unpushed:
            print("  Unpushed assets:")
            for link in unpushed:
                print(f"    - {link.asset_digest} (mode={link.mode})")


REMOTE_DISPATCH: dict[type, Callable] = {
    RemoteListModel: handle_remote_list,
    RemoteAddModel: handle_remote_add,
    RemoteRemoveModel: handle_remote_remove,
    RemoteStatusModel: handle_remote_status,
}


def handle_remote_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = REMOTE_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown remote subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
