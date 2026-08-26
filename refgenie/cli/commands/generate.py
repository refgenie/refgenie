"""The `generate` command group: models and handlers."""

import sys
from collections.abc import Callable
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand

from refgenie.logger import logger


class GenerateSnakefileModel(BaseModel):
    """generate snakefile: generate a Snakemake file."""

    output_path: Path = Field(
        description="Path to save the generated Snakefile.",
        validation_alias=AliasChoices("o", "output-path"),
    )
    snakefile_template_path: Path | None = Field(
        None,
        description="Path to the Snakefile template.",
        validation_alias=AliasChoices("s", "snakefile-template-path"),
    )


class GenerateParser(BaseModel):
    """Intermediate parser for generate subcommands."""

    snakefile: CliSubCommand[GenerateSnakefileModel] = Field(
        description="Generate a Snakemake file."
    )


def handle_generate_snakefile(cmd, refgenie) -> None:
    from refgenie.snakefile.generate import populate_snakefile_template

    populate_snakefile_template(
        refgenie=refgenie,
        snakefile_output_path=cmd.output_path,
        snakefile_template_path=cmd.snakefile_template_path,
    )


GENERATE_DISPATCH: dict[type, Callable] = {
    GenerateSnakefileModel: handle_generate_snakefile,
}


def handle_generate_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = GENERATE_DISPATCH.get(type(leaf))
    if handler is None:
        logger.error(f"Unknown generate subcommand: {type(leaf).__name__}")
        sys.exit(1)
    handler(leaf, refgenie)
