"""The `alias` command group: models and handlers."""

from collections.abc import Callable

from pydantic import AliasChoices, BaseModel, ConfigDict, Field, model_validator
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.commands.framework import CliList
from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.logger import logger


class AliasGetNestedModel(BaseModel):
    """alias get: retrieve aliases or digests (mutually exclusive)."""

    # Without populate_by_name, a field whose validation_alias omits its own name
    # (e.g. genome_digests -> "g"/"genome-digests") cannot be set by field name:
    # pydantic drops the keyword as an unknown extra, so the mutual-exclusion
    # validator never fires for programmatic construction.
    model_config = ConfigDict(populate_by_name=True)

    aliases: CliList | None = Field(
        None,
        description="Aliases to get the digests for.",
        validation_alias=AliasChoices("a", "aliases"),
    )
    genome_digests: CliList | None = Field(
        None,
        description="Genome digests to get the aliases for.",
        validation_alias=AliasChoices("g", "genome-digests"),
    )

    @model_validator(mode="after")
    def check_aliases_digests_exclusion(self):
        if self.aliases is not None and self.genome_digests is not None:
            raise ValueError("--aliases and --genome-digests are mutually exclusive")
        return self


class AliasSetModel(BaseModel):
    """alias set: set genome aliases."""

    aliases: CliList = Field(
        description="Aliases to set.",
        validation_alias=AliasChoices("a", "aliases"),
    )
    digest: str | None = Field(
        None,
        description="Digest to set.",
        validation_alias=AliasChoices("d", "digest"),
    )
    reset: bool = Field(
        False,
        description="Remove all aliases before setting new ones.",
        validation_alias=AliasChoices("r", "reset"),
    )
    force: bool = Field(
        False,
        description="Force action if genome does not exist.",
        validation_alias=AliasChoices("f", "force"),
    )


class AliasRemoveModel(BaseModel):
    """alias remove: remove genome aliases."""

    aliases: CliList = Field(
        description="Aliases to remove.",
        validation_alias=AliasChoices("a", "aliases"),
    )


class AliasParser(BaseModel):
    """Intermediate parser for alias subcommands."""

    get: CliSubCommand[AliasGetNestedModel] = Field(description="Get aliases.")
    set: CliSubCommand[AliasSetModel] = Field(description="Set aliases.")
    remove: CliSubCommand[AliasRemoveModel] = Field(description="Remove aliases.")


def handle_alias_get(cmd, refgenie) -> None:
    rprint(refgenie.alias.table(aliases=cmd.aliases, genome_digests=cmd.genome_digests))


def handle_alias_set(cmd, refgenie) -> None:
    if cmd.digest is not None and not refgenie.genome.exists(cmd.digest):
        if not cmd.force:
            fail(
                f"Genome with digest {cmd.digest} does not exist. "
                "You must initialize it first by building/pulling an asset for that genome.",
                EXIT_NOT_FOUND,
            )
        logger.warning(
            f"--force: genome {cmd.digest} does not exist; creating a dangling alias."
        )
    if cmd.reset:
        for alias in refgenie.alias.get_for_genome(genome_digest=cmd.digest):
            refgenie.alias.remove(alias)
    for alias in cmd.aliases:
        refgenie.set_genome_alias(alias_name=alias, genome_digest=cmd.digest)


def handle_alias_remove(cmd, refgenie) -> None:
    for alias in cmd.aliases:
        refgenie.alias.remove(alias)


ALIAS_DISPATCH: dict[type, Callable] = {
    AliasGetNestedModel: handle_alias_get,
    AliasSetModel: handle_alias_set,
    AliasRemoveModel: handle_alias_remove,
}


def handle_alias_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = ALIAS_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown alias subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
