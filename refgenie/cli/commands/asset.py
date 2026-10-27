"""The `asset` alias command group (provides `asset list`): model and handler."""

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand

from refgenie.cli.commands.listing import handle_list
from refgenie.cli.framework import CliList
from refgenie.cli.errors import fail


class AssetListNestedModel(BaseModel):
    """asset list: list assets (alias of top-level 'list')."""

    # CliList, not list[str]: handle_asset_group forwards this straight to
    # handle_list, so it must accept exactly what ListModel accepts.
    genome: CliList | None = Field(
        default=None,
        description="One or more genome aliases to list assets for.",
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: CliList | None = Field(
        default=None,
        description="Genome digests, used in place of aliases.",
        validation_alias=AliasChoices("genome-digest"),
    )


class AssetParser(BaseModel):
    """Intermediate parser for asset subcommands (alias group)."""

    list: CliSubCommand[AssetListNestedModel] = Field(
        description="List assets (alias of top-level 'list')."
    )


def handle_asset_group(cmd, refgenie) -> None:
    """Top-level 'asset' alias group dispatcher (provides `asset list`)."""
    leaf = get_subcommand(cmd, is_required=True)
    if isinstance(leaf, AssetListNestedModel):
        # Defer to top-level list handler.
        handle_list(leaf, refgenie)
    else:
        fail(f"Unknown asset subcommand: {type(leaf).__name__}")
