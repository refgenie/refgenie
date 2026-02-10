from __future__ import annotations

import argparse
import os
import signal
import sys
from collections.abc import Callable
from logging import getLogger
from time import gmtime, strftime
from typing import Any

import pypiper
from attmap import AttMap
from refgenconf import RefGenConf
from refgenconf.const import (
    BUILD_STATS_DIR,
    CFG_ALIASES_KEY,
    CFG_ASSET_CUSTOM_PROPS_KEY,
    CFG_ASSET_DATE_KEY,
    CFG_ASSET_DESC_KEY,
    CFG_ASSET_PATH_KEY,
    TEMPLATE_TARGET,
)
from refgenconf.exceptions import RefgenconfError
from refgenconf.helpers import block_iter_repr
from ubiquerg import checksum

from .asset_build_packages import *
from .const import *
from .helpers import _writeable, make_sure_path_exists

_LOGGER = getLogger(PKG_NAME)


def _handle_sigint(gat: list[str]) -> Callable[[int, Any], None]:
    """Create a SIGINT handler that logs the interruption and exits.

    Args:
        gat: A list of [genome, asset, tag] for message generation.

    Returns:
        The SIGINT handling function.
    """

    def handle(sig, frame):
        _LOGGER.warning("\nThe build was interrupted: {}/{}:{}".format(*gat))
        sys.exit(0)

    return handle


def _seek(
    rgc: RefGenConf,
    genome_name: str,
    asset_name: str,
    tag_name: str | None = None,
    seek_key: str | None = None,
    enclosing_dir: bool = False,
    strict_exists: bool = False,
) -> str:
    """Seek with strict file existence check."""
    return rgc.seek_src(
        genome_name=genome_name,
        asset_name=asset_name,
        tag_name=tag_name,
        seek_key=seek_key,
        enclosing_dir=enclosing_dir,
        strict_exists=strict_exists,
    )


def test_recipe(
    rgc: RefGenConf,
    recipe: Any,
    asset_class: Any,
    test_inputs: dict[str, Any],
    alias: str,
    args: argparse.Namespace,
) -> tuple[bool, str | None]:
    """Test a recipe by building it and verifying outputs."""
    _LOGGER.info(f"Starting {recipe.name} recipe test")
    _LOGGER.info(f"Test inputs:\n{test_inputs}")
    asset = recipe.name
    genome = recipe.test["genome"]
    genome_outfolder = os.path.join(rgc.data_dir, genome)
    build_namespaces = AttMap(
        {
            "genome": genome,
            "asset": asset,
            "params": test_inputs["params"],
            "files": test_inputs["files"],
            "custom_properties": recipe.resolve_custom_properties(
                use_docker=args.docker
            ),
            "assets": test_inputs["assets"],
        }
    )

    tag = recipe.resolve_default_tag(namespaces=build_namespaces)
    build_namespaces.update({"tag": tag})
    build_namespaces.update(
        {"asset_outfolder": os.path.join(genome_outfolder, asset, tag)}
    )
    # rgc.make_writable(
    #     filepath=os.path.join(rgc.data_dir, "_recipe_test", "temp_genome_config.yaml")
    # )

    build_stats_dir = os.path.abspath(
        os.path.join(genome_outfolder, asset, tag, BUILD_STATS_DIR)
    )
    if not os.path.exists(build_stats_dir):
        os.makedirs(build_stats_dir, exist_ok=True)
    _LOGGER.info(
        f"Saving outputs to:{block_iter_repr([f'content: {genome_outfolder}', f'logs: {build_stats_dir}'])}"
    )

    if not _writeable(genome_outfolder):
        raise OSError(
            f"Insufficient permissions to write to output folder: {genome_outfolder}"
        )

    pipeline_name = (
        f"{PKG_NAME}_build_{build_namespaces['asset']}_{build_namespaces['tag']}"
    )
    pm = pypiper.PipelineManager(
        name=pipeline_name,
        outfolder=build_stats_dir,
        args=args,
    )
    _LOGGER.debug("Recipe: " + str(recipe))

    # create a bundle list to simplify calls below
    gat = [genome, asset, build_namespaces["tag"]]

    # populate recipe commands with build_namespaces
    command_list_populated = recipe.populate_command_templates(
        namespaces=build_namespaces
    )

    # create output directory if it doesn't exist
    make_sure_path_exists(build_namespaces["asset_outfolder"])

    target = os.path.join(build_stats_dir, TEMPLATE_TARGET.format(genome, asset, tag))
    # add target command
    command_list_populated.append(f"touch {target}")
    try:
        # run build command
        signal.signal(signal.SIGINT, _handle_sigint(gat))
        return_code = pm.run(
            command_list_populated,
            target,
            container=pm.container,
            default_return_code=None,
        )
    except pypiper.exceptions.SubprocessError:
        return False, f"Asset '{genome}/{asset}:{tag}' build failed"
    else:
        # add updates to config file
        with rgc as r:
            if asset == "fasta":
                r.update_genomes(
                    genome, data={CFG_ALIASES_KEY: [alias]}, force_digest=genome
                )
            r.update_assets(
                *gat[0:2],
                data={CFG_ASSET_DESC_KEY: recipe.description},
                force_digest=genome,
            )
            try:
                r.set_asset_class(genome, asset, recipe.output_class.name)
            except RefgenconfError:
                _LOGGER.error(
                    "You can't mix assets of different classes within a single asset namespace"
                )
                raise
            r.update_tags(
                *gat,
                force_digest=genome,
                data={
                    CFG_ASSET_PATH_KEY: asset,
                    CFG_ASSET_CUSTOM_PROPS_KEY: build_namespaces["custom_properties"],
                    CFG_ASSET_DATE_KEY: strftime("%Y-%m-%d_%H:%M", gmtime()),
                },
            )
            r.update_seek_keys(
                *gat,
                force_digest=genome,
                keys={
                    k: v.format(**build_namespaces)
                    for k, v in recipe.output_class.seek_keys.items()
                },
            )
            r.set_default_pointer(*gat, force_digest=genome)
            r.set_genome_alias(genome=alias, digest=genome)
    pm.stop_pipeline()
    test_outputs = recipe.get_test_outputs()
    # run tests
    _LOGGER.info(f"Testing '{recipe.name}' recipe outputs")
    for seek_key, output_values_tests in test_outputs.items():
        if seek_key not in asset_class.seek_keys:
            raise ValueError(
                f"Seek key '{seek_key}' cannot be tested. Seek keys defined in '{recipe.output_class.name}' "
                f"asset class: {', '.join(asset_class.seek_keys.keys())}"
            )
        value_to_test = _seek(rgc, genome, asset, tag, seek_key)
        _LOGGER.info(f"Testing result: {value_to_test} ({seek_key})")
        if "md5sum" in output_values_tests:
            if not os.path.exists(value_to_test):
                raise OSError(
                    f"Nothing to run checksum on, '{value_to_test} does not exist"
                )
            actual_checksum = checksum(value_to_test)
            if output_values_tests["md5sum"] != actual_checksum:
                return (
                    False,
                    f"Checksum mismatch for '{seek_key}' ({output_values_tests['md5sum']} != {actual_checksum})",
                )
        if "value" in output_values_tests:
            if output_values_tests["value"] != value_to_test:
                return (
                    False,
                    f"Value mismatch for '{seek_key}' ({output_values_tests['value']} != {value_to_test})",
                )
    return True, None
