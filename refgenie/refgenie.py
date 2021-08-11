import csv
import json
import os
import re
import signal
import sys
from glob import glob
from logging import getLogger
from time import gmtime, strftime

import attmap
import pypiper
from refgenconf import RefGenConf, get_dir_digest, recipe_factory
from refgenconf.const import CFG_ASSET_TAGS_KEY, CFG_CHECKSUM_KEY, TEMPLATE_RECIPE_JSON
from refgenconf.exceptions import (
    MissingAssetError,
    MissingGenomeError,
    MissingSeekKeyError,
    MissingTagError,
    RefgenconfError,
)
from refgenconf.helpers import block_iter_repr
from rich.console import Console
from rich.progress import track
from ubiquerg import is_url
from ubiquerg import parse_registry_path as prp
from ubiquerg.files import checksum
from ubiquerg.system import is_writable
from yacman import UndefinedAliasError

from .asset_build_packages import *
from .const import *
from .helpers import (
    _parse_user_kw_input,
    _raise_missing_recipe_error,
    _skip_lock,
    _writeable,
    make_sure_path_exists,
)

_LOGGER = getLogger(PKG_NAME)


def parse_registry_path(path):
    return prp(
        path,
        defaults=[
            ("protocol", None),
            ("genome", None),
            ("asset", None),
            ("seek_key", None),
            ("tag", None),
        ],
    )


def get_build_namespaces(
    genome,
    asset,
    tag,
    outfolder,
    specified_files=None,
    specified_params=None,
    **kwargs,
):
    """
    Gives a dict with variables used to populate an build commands with.
    """
    build_namespaces = attmap.AttMap(
        {
            "genome": genome,
            "asset": asset,
            "tag": tag,
            "asset_outfolder": os.path.join(outfolder, asset, tag),
            "params": specified_params,
            "files": specified_files,
        }
    )
    build_namespaces.update(**kwargs)
    return build_namespaces


def refgenie_initg(rgc, genome, content_checksums):
    """
    Initializing a genome means adding `collection_checksum` attributes in the
    genome config file. This should perhaps be a function in refgenconf, but not
    a CLI-hook. Also adds `content_checksums` tsv file (should be a recipe cmd?).

    This function updates the provided RefGenConf object with the
    genome(collection)-level checksum and saves the individual checksums to a
    TSV file in the fasta asset directory.

    :param refgenconf.RefGenConf rgc: genome configuration object
    :param str genome: name of the genome
    :param dict content_checksums: checksums of individual content_checksums, e.g. chromosomes
    """
    genome_dir = os.path.join(rgc.data_dir, genome)
    if is_writable(genome_dir):
        output_file = os.path.join(genome_dir, "{}_sequence_digests.tsv".format(genome))
        with open(output_file, "w") as contents_file:
            wr = csv.writer(contents_file, delimiter="\t")
            for key, val in content_checksums.items():
                wr.writerow([key, val])
        _LOGGER.debug("sequence digests saved to: {}".format(output_file))
    else:
        _LOGGER.warning(
            "Could not save the genome sequence digests. '{}' is not writable".format(
                genome_dir
            )
        )


def refgenie_build_reduce(gencfg, preserve_map_configs=False):
    """
    Asset building process may be split into two tasks: building assets (_Map_ procedure)
    and gathering asset metadata (_Reduce_ procedure).

    This function performs the _Reduce_ procedure:
    finds the genome configuration files produced in the _Map_ step,
    updates the main genome configuration file with their contents and removes them.

    :param str gencfg: an absolute path to the genome configuration file
    :param bool preserve_map_configs: a boolean indicating whether the map configs should be preserved,
        by default they are removed once the contents are integrated into the master genome config.
    :return bool: a boolean indicating whether the master config has been successfully updated
        or None in case there were no map configs found.
    """

    def _map_cfg_match_pattern(data_dir, match_all_str):
        """
        Create a path to the map genome config witb a provided 'match all' character,
        which needs to be different depending on the matchig scenario.

        :param str data_dir: an absolute path to the data directory
        :param str match_all_str: match all character to use
        """
        return os.path.join(
            data_dir,
            *([match_all_str] * 3),
            BUILD_STATS_DIR,
            BUILD_MAP_CFG,
        )

    _LOGGER.info("Running the reduce procedure. No assets will be built.")
    rgc_master = RefGenConf(filepath=gencfg, writable=True)
    regex_pattern = _map_cfg_match_pattern(rgc_master.data_dir, "(\S+)")
    glob_pattern = _map_cfg_match_pattern(rgc_master.data_dir, "*")
    rgc_map_filepaths = glob(glob_pattern, recursive=True)
    if len(rgc_map_filepaths) == 0:
        _LOGGER.info(f"No map configs to reduce")
        return None
    _LOGGER.debug(f"Map configs to reduce: {block_iter_repr(rgc_map_filepaths)}")
    matched_gats = []
    for rgc_map_filepath in track(
        rgc_map_filepaths,
        description=f"Reducing {len(rgc_map_filepaths)} configs",
    ):
        matched_genome, matched_asset, matched_tag = re.match(
            pattern=regex_pattern, string=rgc_map_filepath
        ).groups()
        matched_gat = f"{matched_genome}/{matched_asset}:{matched_tag}"
        map_rgc = RefGenConf(filepath=rgc_map_filepath, writable=False)
        if CFG_GENOMES_KEY not in map_rgc:
            _LOGGER.warning(
                f"'{rgc_map_filepath}' is missing '{CFG_GENOMES_KEY}' key, skipping"
            )
            continue
        # this should be a one element list
        genome_digests = map_rgc[CFG_GENOMES_KEY].keys()
        if len(genome_digests) > 1:
            _LOGGER.warning(
                f"There are {len(genome_digests)} genomes in the map build config while 1 expected, skipping"
            )
            continue
        genome_digest = genome_digests[0]
        alias = map_rgc.get_genome_alias(digest=genome_digest)
        if genome_digest != matched_genome:
            raise Exception(
                f"Genome directory name does not match genome in the map config: {matched_genome} != {genome_digest}"
            )
        asset_data = tag_data = map_rgc[CFG_GENOMES_KEY][matched_genome][
            CFG_ASSETS_KEY
        ][matched_asset]
        tag_data = asset_data[CFG_ASSET_TAGS_KEY][matched_tag]
        default_tag_in_map = asset_data[CFG_ASSET_DEFAULT_TAG_KEY]
        try:
            alias_master = rgc_master.get_genome_alias(digest=genome_digest)
            assert alias == alias_master
        except (UndefinedAliasError, AssertionError):
            # no need to put this in context manager
            # it is already used in the method
            rgc_master.set_genome_alias(
                genome=alias, digest=genome_digest, create_genome=True
            )
        with rgc_master as r:
            if CFG_ASSET_PARENTS_KEY in tag_data:
                for parent in tag_data[CFG_ASSET_PARENTS_KEY]:
                    parsed_parent = parse_registry_path(parent)
                    r.update_relatives_assets(
                        genome=parsed_parent["genome"],
                        asset=parsed_parent["asset"],
                        tag=parsed_parent["tag"],
                        data=[matched_gat],
                        children=True,
                    )

            if CFG_ASSET_CHILDREN_KEY in tag_data:
                for child in tag_data[CFG_ASSET_CHILDREN_KEY]:
                    parsed_child = parse_registry_path(child)
                    r.update_relatives_assets(
                        genome=parsed_child["genome"],
                        asset=parsed_child["asset"],
                        tag=parsed_child["tag"],
                        data=[matched_gat],
                        children=False,
                    )
            r.update_tags(
                genome=matched_genome,
                asset=matched_asset,
                tag=matched_tag,
                data=tag_data,
                force_digest=genome_digest,
            )
            # set a default tag in the master config to the one built in map mode,
            # this will not overwrite an existing tag though
            r.set_default_pointer(
                genome=matched_genome,
                asset=matched_asset,
                tag=default_tag_in_map,
            )
            r.set_asset_class(
                genome=matched_genome,
                asset=matched_asset,
                asset_class=asset_data[CFG_ASSET_CLASS_KEY],
            )
        matched_gats.append(matched_gat)
        if not preserve_map_configs:
            os.remove(rgc_map_filepath)
    _LOGGER.info(f"Added entries for: {block_iter_repr(matched_gats)}")
    return True


def refgenie_build(gencfg, genome, asset_list, recipe_source, args, pipeline_kwargs):
    """
    Runs the refgenie build recipe.

    :param str gencfg: path to the genome configuration file
    :param str genome: the genome to build
    :param list asset_list: a list of assets to build
    :param str recipe_source: the name of the recipe to use to build the assets
    :param argparse.Namespace args: parsed command-line options/arguments
    """

    rgc = RefGenConf(
        filepath=gencfg,
        writable=False,
        skip_read_lock=_skip_lock(args.skip_read_lock, gencfg),
    )
    specified_files = _parse_user_kw_input(args.files)
    specified_params = _parse_user_kw_input(args.params)

    def _build_asset(
        build_namespaces,
        recipe,
        alias,
        pipeline_kwargs,
    ):
        """
        Builds assets with pypiper and updates a genome config file.

        This function actually runs the build commands in a given build package,
        and then update the refgenie config file.

        :param attmap.AttMap build_namespaces: a mapping of namespaces to populate the templates with
        :param str alias: the genome alias to use
        :param dict pipeline_kwargs: the kwargs to pass to the pypiper pipeline
        :param refgenconf.Recipe recipe: A recipe object specifying lists
            of required input_assets, commands to run, and outputs to register as
            assets.
        """
        if args.map:
            # Performing a build map step.
            # The reduce step will need to be performed to get the built
            # asset metadata to the master config file
            genome_alias = rgc.get_genome_alias(digest=genome)
            # create an empty config file in the genome directory
            _LOGGER.info(f"Using new map genome config: {locked_map_gencfg}")
            make_sure_path_exists(os.path.dirname(locked_map_gencfg))
            open(locked_map_gencfg, "a").close()
            # initialize a new RefGenConf.
            # Use the master location for data storage,
            # but change path to the in asset dir location
            rgc_map = RefGenConf(
                entries={"genome_folder": rgc.genome_folder},
                filepath=locked_map_gencfg,
            )
            # set the alias first (if available), based on the master file

            rgc_map.set_genome_alias(
                digest=genome,
                genome=genome_alias,
                create_genome=True,
            )

            # copy the genome of interest section to the new RefGenConf,
            # so that possible dependancies can be satisfied
            rgc_map.update_genomes(
                genome=genome_alias,
                data=rgc[CFG_GENOMES_KEY][genome],
            )

        else:
            rgc_map = rgc

        _LOGGER.info(
            f"Saving outputs to:{block_iter_repr([f'content: {genome_outfolder}', f'logs: {build_stats_dir}'])}"
        )

        if not _writeable(genome_outfolder):
            _LOGGER.error(
                f"Insufficient permissions to write to output folder: {genome_outfolder}"
            )
            return False, rgc_map

        pipeline_name = (
            pipeline_kwargs["pipeline_name"]
            if "pipeline_name" in pipeline_kwargs
            else None
            or f"{PKG_NAME}_build_{build_namespaces['asset']}_{build_namespaces['tag']}"
        )
        pm = pypiper.PipelineManager(
            name=pipeline_name,
            outfolder=build_stats_dir,
            args=args,
            **pipeline_kwargs,
        )
        if args.docker:
            # Set up some docker stuff
            volumes = args.volumes or []
            volumes.append(genome_outfolder)
            pm.get_container(recipe.container, volumes)
        _LOGGER.debug("Recipe: " + str(recipe))

        inputs = {
            "files": specified_files,
            "params": specified_params,
            "assets": input_assets_dict,
        }
        # save inputs to json file
        inputs_file = rgc.get_recipe_inputs_path(genome, asset, tag)
        with open(inputs_file, "w") as f:
            json.dump(inputs, f, indent=4)
        _LOGGER.debug(f"Using inputs: {inputs}\nInputs saved to: {inputs_file}")

        # create a bundle list to simplify calls below
        gat = [genome, asset, build_namespaces["tag"]]

        # populate recipe commands with build_namespaces
        command_list_populated = recipe.populate_command_templates(
            namespaces=build_namespaces
        )

        # create output directory if it doesn't exist
        make_sure_path_exists(build_namespaces["asset_outfolder"])

        target = os.path.join(
            build_stats_dir, TEMPLATE_TARGET.format(genome, asset, tag)
        )
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
            _LOGGER.error(f"Asset '{genome}/{asset}:{tag}' build failed")
            return False, rgc_map
        else:
            if return_code is None:
                # if no commands were run, stop the pipeline and return
                pm.stop_pipeline()
                return None, rgc_map
            # save build recipe to the JSON-formatted file
            recipe_file_path = os.path.join(
                build_stats_dir, TEMPLATE_RECIPE_JSON.format(asset, tag)
            )
            recipe.to_json(filepath=recipe_file_path)
            # since the assets are always built to a standard dir structure, we
            # can just stitch a path together for asset digest calculation
            asset_dir = os.path.join(rgc_map.data_dir, *gat)
            if not os.path.exists(asset_dir):
                raise OSError(
                    f"Could not compute asset digest. Path does not exist: {asset_dir}"
                )
            digest = get_dir_digest(
                path=asset_dir, exclude_files=recipe.checksum_exclude_list
            )
            _LOGGER.info(f"Asset digest: {digest}")
            # add updates to config file
            with rgc_map as r:
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
                        f"You can't mix assets of different classes within a single asset namespace"
                    )
                    raise
                r.update_tags(
                    *gat,
                    force_digest=genome,
                    data={
                        CFG_ASSET_PATH_KEY: asset,
                        CFG_ASSET_CHECKSUM_KEY: digest,
                        CFG_ASSET_CUSTOM_PROPS_KEY: build_namespaces[
                            "custom_properties"
                        ],
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
        pm.stop_pipeline()
        return True, rgc_map

    for single_asset in asset_list:
        asset = single_asset["asset"]
        recipe_source = recipe_source or asset
        recipe = rgc.get_recipe(recipe_name=recipe_source)
        # handle user-requested parents for the required assets
        assets = {}
        input_assets_dict = {}
        parent_assets = []
        specified_asset_keys, specified_assets = None, None
        if args.assets is not None:
            parsed_parents_input = _parse_user_kw_input(args.assets)
            specified_asset_keys = list(parsed_parents_input.keys())
            specified_assets = list(parsed_parents_input.values())
            _LOGGER.info(
                f"Custom assets requested: {block_iter_repr(args.assets, flatten=True)}"
            )
        if not specified_asset_keys and isinstance(args.assets, list):
            _LOGGER.warning(
                "Specified parent assets format is invalid. Using defaults."
            )
        for recipe_name, req_asset_map in recipe.required_assets.items():
            # for each req asset see if non-default parents were requested
            req_class_name = req_asset_map["asset_class"]
            if (
                specified_asset_keys is not None
                and req_class_name in specified_asset_keys
            ):
                spec_idx = specified_asset_keys.index(req_class_name)
                parent_data = parse_registry_path(specified_assets[spec_idx])
                g, a, t, s = (
                    parent_data["genome"],
                    parent_data["asset"],
                    parent_data["tag"]
                    or rgc.get_default_tag(genome, parent_data["asset"]),
                    specified_asset_keys[spec_idx],
                )
            else:  # if no custom parents requested for the req asset, use default one
                default = parse_registry_path(req_asset_map[DEFAULT])
                g, a, t, s = (
                    genome,
                    default["asset"],
                    default["tag"] or rgc.get_default_tag(genome, default["asset"]),
                    default["seek_key"],
                )
            effective_asset_class = rgc.get_assets_asset_class(g, a)
            if req_class_name != effective_asset_class:
                raise RefgenconfError(
                    f"Class of the input asset ({g}/{a}:{t}) does not match the"
                    f"'{recipe_name}' recipe input class requirement: "
                    f"{effective_asset_class} != {req_class_name}"
                )
            try:
                parent_assets.append(
                    f"{rgc.get_genome_alias_digest(g, fallback=True)}/{a}:{t}"
                )
            except UndefinedAliasError as e:
                _LOGGER.warning(f"'{g}' namespace has not been initialized yet")
                if args.pull_parents:
                    _LOGGER.info(f"Pulling missing parent: {g}/{a}:{t}")
                    ret = rgc.pull(genome=g, asset=a, tag=t)
                    if ret is None or not all(ret):
                        _LOGGER.info(
                            f"Missing parent asset pull requested, but failed: {g}/{a}:{t}. "
                            f"Reason: {str(e)}"
                        )
                        return False
                else:
                    raise
            try:
                assets[req_class_name] = _seek(rgc, g, a, t, s)
                input_assets_dict[req_class_name] = parent_assets[0]
            except (
                MissingAssetError,
                MissingGenomeError,
                MissingTagError,
                MissingSeekKeyError,
            ) as e:
                if args.pull_parents:
                    _LOGGER.info(f"Pulling missing parent: {g}/{a}:{t}")
                    ret = rgc.pull(genome=g, asset=a, tag=t)
                    if ret is None or not all(ret):
                        _LOGGER.info(
                            f"Missing parent asset pull requested, but failed: {g}/{a}:{t}. "
                            f"Reason: {str(e)}"
                        )
                        return False
                else:
                    raise

        for required_file_name, required_file in recipe.required_files.items():
            if (
                specified_files is None
                or required_file_name not in specified_files.keys()
            ):
                raise ValueError(
                    "Path to the '{x}' input ({desc}) is required, but not provided. "
                    "Specify it with: --files {x}=/path/to/{x}_file".format(
                        x=required_file_name, desc=required_file[DESC]
                    )
                )
        for required_param_name, required_param in recipe.required_params.items():
            if specified_params is None:
                specified_params = {}
            if required_param_name not in specified_params.keys():
                if required_param[DEFAULT] is None:
                    raise ValueError(
                        "Value for the parameter '{x}' ({desc}) is required, but not provided. "
                        "Specify it with: --params {x}=value".format(
                            x=required_param_name, desc=required_param[DESC]
                        )
                    )
                else:
                    specified_params.update(
                        {required_param_name: required_param[DEFAULT]}
                    )
        ori_genome = genome
        if recipe.output_class.name == "fasta":
            if genome in rgc.genomes_list() and "fasta" in rgc.list_assets_by_genome(
                genome
            ):
                pretag = rgc.get_default_tag(genome, "fasta")
                _LOGGER.warning(
                    "'{g}' genome is already initialized with other fasta asset ({g}/{a}:{t})".format(
                        g=genome, a=asset, t=pretag
                    )
                )
                genome = rgc.get_genome_alias_digest(alias=genome, fallback=True)
            else:
                # if the recipe is "fasta" we first initialiaze the genome, based on the provided path to the input FASTA file
                genome, _ = rgc.initialize_genome(
                    fasta_path=specified_files["fasta"],
                    alias=ori_genome,
                    skip_alias_write=True,
                )
        else:
            try:
                genome = rgc.get_genome_alias_digest(genome, fallback=True)
            except UndefinedAliasError:
                _LOGGER.error(
                    f"Genome '{genome}' has not been initialized yet; "
                    "no key found for this alias"
                )
                return False

        genome_outfolder = os.path.join(rgc.data_dir, genome)
        # collect variables required to populate the command templates
        build_namespaces = attmap.AttMap(
            {
                "genome": genome,
                "asset": asset,
                "params": specified_params,
                "files": specified_files,
                "custom_properties": recipe.resolved_custom_properties,
                "assets": assets,
            }
        )
        # determine tag; it's either forced or defined in the recipe or default
        tag = (
            single_asset["tag"]
            or recipe.resolve_default_tag(namespaces=build_namespaces)
            or rgc.get_default_tag(genome, asset, use_existing=False)
        )
        if any([c in tag for c in TAG_NAME_BANNED_CHARS]):
            raise ValueError(
                f"The tag name can't consist of characters: {TAG_NAME_BANNED_CHARS}"
            )
        build_namespaces["tag"] = tag
        build_namespaces.update(
            {"asset_outfolder": os.path.join(genome_outfolder, asset, tag)}
        )
        build_stats_dir = os.path.abspath(
            os.path.join(genome_outfolder, asset, tag, BUILD_STATS_DIR)
        )
        locked_map_gencfg = os.path.join(build_stats_dir, LOCKED_BUILD_MAP_CFG)
        map_gencfg = os.path.join(build_stats_dir, BUILD_MAP_CFG)

        _LOGGER.info(
            f"Building '{genome}/{asset}:{tag}' using '{recipe_source}' recipe"
        )
        is_built, rgc_map = _build_asset(
            build_namespaces,
            recipe,
            ori_genome,
            pipeline_kwargs,
        )
        if is_built is None:
            _LOGGER.info("No commands were executed")
            return True
        if not is_built:
            log_path = os.path.abspath(
                os.path.join(
                    genome_outfolder,
                    asset,
                    tag,
                    BUILD_STATS_DIR,
                    ORI_LOG_NAME_REGEX,
                )
            )
            _LOGGER.info(
                f"'{genome}/{asset}:{tag}' was not added to the config, but directory has been left in place. "
                f"See the log file for details: {log_path}"
            )
            return False
        _LOGGER.info(f"Finished building '{asset}' asset")
        with rgc_map as r:
            # update asset relationships
            r.update_relatives_assets(genome, asset, tag, parent_assets)  # adds parents
            for i in parent_assets:
                parsed_parent = parse_registry_path(i)
                # adds child (currently built asset) to the parent
                r.update_relatives_assets(
                    parsed_parent["genome"],
                    parsed_parent["asset"],
                    parsed_parent["tag"],
                    [f"{genome}/{asset}:{tag}"],
                    True,
                )
            if args.genome_description is not None:
                _LOGGER.debug(
                    f"Adding genome ({genome}) description: '{args.genome_description}'"
                )
                r.update_genomes(genome, {CFG_GENOME_DESC_KEY: args.genome_description})
            if args.tag_description is not None:
                _LOGGER.debug(
                    f"Adding tag ({genome}/{asset}:{tag}) description: '{args.tag_description}'"
                )
                r.update_tags(
                    genome,
                    asset,
                    tag,
                    {CFG_TAG_DESC_KEY: args.tag_description},
                )
        rgc_map._symlink_alias(genome, asset, tag)
        if args.map:
            # move the contents of the locked map config to a map config,
            # which is discoverable by the reduce step
            os.rename(locked_map_gencfg, map_gencfg)
            _LOGGER.info(
                f"Asset metadata saved in '{map_gencfg}'. "
                f"To make the asset accessible globally run 'refgenie build --reduce'"
            )

        return True


def _handle_sigint(gat):
    """
    SIGINT handler, unlocks the config file and exists the program

    :param list gat: a list of genome, asset and tag. Used for a message generation.
    :return function: the SIGINT handling function
    """

    def handle(sig, frame):
        _LOGGER.warning("\nThe build was interrupted: {}/{}:{}".format(*gat))
        sys.exit(0)

    return handle


def _seek(
    rgc, genome_name, asset_name, tag_name=None, seek_key=None, enclosing_dir=False
):
    """
    Strict seek. Most use cases in this package require file existence
     check in seek. This function makes it easier
    """
    return rgc.seek_src(
        genome_name=genome_name,
        asset_name=asset_name,
        tag_name=tag_name,
        seek_key=seek_key,
        enclosing_dir=enclosing_dir,
        strict_exists=True,
    )
