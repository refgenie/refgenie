
import glob
import os
import yacman
from logging import getLogger
from refgenconf import (
    DownloadJsonError,
    MissingAssetError,
    MissingGenomeError,
    RefGenConf,
)

from .helpers import (
    _parse_user_build_input,
    _raise_missing_recipe_error,
    _skip_lock,
    _writeable,
)

from .const import *
from .refgenie import get_asset_vars


_LOGGER = getLogger(PKG_NAME)


def get_flexible_recipes():
    asset_class_files = []
    for file in glob.glob("asset_classes/*.yaml"):
        asset_class_files.append(file)

    asset_classes = {}
    for ac_path in asset_class_files:
        aclass = yacman.YacAttMap(filepath=ac_path)
        asset_classes[aclass["class_name"]] = aclass

    recipe_files = []
    for file in glob.glob("recipes/*.yaml"):
        recipe_files.append(file)

    recipes = {}
    for r_path in recipe_files:
        recipe = yacman.YacAttMap(filepath=r_path)
        recipe["output"] = asset_classes[recipe["output_class"]]
        recipes[recipe["recipe_name"]] = recipe

    return {"recipes" : recipes, "asset_classes": asset_classes}


def refgenie_build_flex(gencfg, genome, asset_list, recipe_name, args):

    rgc = RefGenConf(
        filepath=gencfg,
        writable=False,
        skip_read_lock=_skip_lock(args.skip_read_lock, gencfg),
    )

    if not hasattr(args, "outfolder") or not args.outfolder:
        # Default to genome_folder
        _LOGGER.debug("No outfolder provided, using genome config.")
        args.outfolder = rgc.data_dir


    specified_args = _parse_user_build_input(args.files)
    specified_params = _parse_user_build_input(args.params)

    for a in asset_list:
        asset_key = a["asset"]
        asset_tag = a["tag"] or rgc.get_default_tag(
            genome, a["asset"], use_existing=False
        )

    genome_outfolder = os.path.join(args.outfolder, genome)

    asset_vars = get_asset_vars(
        genome,
        asset_key,
        asset_tag,
        genome_outfolder,
        specified_args,
        specified_params
    )

    flex_recipes = get_flexible_recipes()

    print(flex_recipes)
    print(asset_vars)

    return None
