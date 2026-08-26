import os
from pathlib import Path

from jinja2 import Template
from pydantic import BaseModel

from refgenie import Refgenie
from refgenie.logger import logger

SNAKEFILE_TEMPLATE_THREADS = int(os.environ.get("SNAKEFILE_TEMPLATE_THREADS", 4))
SNAKEFILE_TEMPLATE_PATH = Path(__file__).parent / "template" / "template.smk"


class AssetBuildRuleSpec(BaseModel):
    asset_name: str
    input_asset_names: list[str]
    input_files: list[str]
    input_params: dict
    needs_genome_init: bool = True


def get_snakefile_context(
    refgenie: Refgenie,
) -> dict[str, list[str] | list[AssetBuildRuleSpec]]:
    """
    Get the context for the Snakemake template.

    Each asset build rule includes only the recipe's actual input_assets (no artificial
    fasta injection). All rules depend on a genome_init rule instead, which ensures the
    genome is initialized before any asset build.

    Args:
        refgenie: The Refgenie instance.

    Returns:
        A dictionary with context data for the Snakemake template.
    """
    # One rule per recipe NAME, not per (name, version): `list_all()` returns every
    # stored version, and snakemake aborts on duplicate rule names; since recipes are
    # immutable, any edit adds a new version. get() with no version resolves the
    # latest. dict.fromkeys preserves order so generated Snakefiles stay diffable.
    recipe_names = list(dict.fromkeys(recipe.name for recipe in refgenie.recipe.list_all()))
    asset_build_rule_specs = []

    for recipe_name in recipe_names:
        recipe = refgenie.recipe.get(recipe_name)
        input_assets = recipe.input_assets or {}
        input_asset_names = list(input_assets.keys())
        input_files = recipe.input_files or {}
        input_params = {"threads": {"default": SNAKEFILE_TEMPLATE_THREADS}}
        input_params.update(recipe.input_params or {})

        asset_build_rule_specs.append(
            AssetBuildRuleSpec(
                asset_name=recipe_name,
                input_asset_names=input_asset_names,
                input_files=list(input_files.keys()),
                input_params={k: v["default"] for k, v in input_params.items()},
                needs_genome_init=True,
            )
        )

    return {
        "asset_names": recipe_names,
        "asset_build_rule_specs": asset_build_rule_specs,
        "has_genome_init_rule": True,
    }


def populate_snakefile_template(
    refgenie: Refgenie,
    snakefile_output_path: Path,
    snakefile_template_path: Path | None = None,
):
    """
    Populate the Snakemake template with the context data and save it to a file.

    Args:
        refgenie: The Refgenie instance.
        snakefile_template_path: The path to the Snakemake template file.
        snakefile_output_path: The path to save the generated Snakemake file.
    """
    snakefile_template_path = snakefile_template_path or SNAKEFILE_TEMPLATE_PATH
    render_context = get_snakefile_context(refgenie)
    logger.info("Rendering the Snakemake template")
    logger.debug(render_context)
    template = Template(snakefile_template_path.read_text())
    logger.info(f"Saving the populated Snakemake template to {snakefile_output_path.resolve()}")
    snakefile_output_path.write_text(template.render(**render_context))

