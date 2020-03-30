#!/usr/bin/env python

from argparse import SUPPRESS
from collections import OrderedDict
from shutil import rmtree
from re import sub
from requests import ConnectionError
import os
import sys
import csv
import signal
import json

import pyfaidx

from ._version import __version__
from .exceptions import MissingGenomeConfigError, MissingFolderError
from .asset_build_packages import *
from .const import *

import logmuse
import pypiper
import refgenconf
from refgenconf import RefGenConf, MissingAssetError, MissingGenomeError, MissingRecipeError, DownloadJsonError
from ubiquerg import is_url, query_yes_no, parse_registry_path as prp, VersionInHelpParser, is_command_callable
from ubiquerg.system import is_writable
import yacman

# from refget import fasta_checksum
from .refget import fasta_checksum

_LOGGER = None


def build_argparser():
    """
    Builds argument parser.

    :return argparse.ArgumentParser
    """

    banner = "%(prog)s - reference genome asset manager"
    additional_description = "\nhttps://refgenie.databio.org"

    parser = VersionInHelpParser(
        prog="refgenie",
        version=__version__,
        description=banner,
        epilog=additional_description)

    subparsers = parser.add_subparsers(dest="command")

    def add_subparser(cmd, description):
        return subparsers.add_parser(
            cmd, description=description, help=description)

    sps = {}
    for cmd, desc in SUBPARSER_MESSAGES.items():
        sps[cmd] = add_subparser(cmd, desc)
        # It's required for init
        sps[cmd].add_argument(
            '-c', '--genome-config', required=(cmd == INIT_CMD), dest="genome_config",
            help="Path to local genome configuration file. Optional if {} environment variable is set."
                .format(", ".join(refgenconf.CFG_ENV_VARS)))

    sps[INIT_CMD].add_argument('-s', '--genome-server', nargs='+', default=DEFAULT_SERVER,
                               help="URL(s) to use for the {} attribute in config file. Default: {}."
                               .format(DEFAULT_SERVER, CFG_SERVERS_KEY))
    sps[BUILD_CMD] = pypiper.add_pypiper_args(
        sps[BUILD_CMD], groups=None, args=["recover", "config", "new-start"])

    # Add any arguments specific to subcommands.

    sps[BUILD_CMD].add_argument(
        '--tag-description', required=False, default=None, type=str,
        help="Add tag level description (e.g. built with version 0.3.2).")

    sps[BUILD_CMD].add_argument(
        '--genome-description', required=False, default=None, type=str,
        help="Add genome level description (e.g. The mouse mitochondrial genome, released in Dec 2013).")

    sps[BUILD_CMD].add_argument(
        "-d", "--docker", action="store_true", help="Run all commands in the refgenie docker container.")

    sps[BUILD_CMD].add_argument(
        '--assets', nargs="+", action='append', required=False, default=None,
        help='Override the default genome, asset and tag of the parents'
             ' (e.g. fasta=hg38/fasta:default gtf=mm10/gencode_gtf:default).')

    sps[BUILD_CMD].add_argument(
        '--files', nargs="+", action='append', required=False, default=None,
        help='Provide paths to the required files (e.g. fasta=/path/to/file.fa.gz).')

    sps[BUILD_CMD].add_argument(
        '--params', nargs="+", action='append', required=False, default=None,
        help='Provide required parameter values (e.g. param1=value1).')

    sps[BUILD_CMD].add_argument(
        '-v', '--volumes', nargs="+", required=False, default=None,
        help='If using docker, also mount these folders as volumes.')

    sps[BUILD_CMD].add_argument(
        '-o', '--outfolder', dest='outfolder', required=False, default=None,
        help='Override the default path to genomes folder, which is the '
             'genome_folder attribute in the genome configuration file.')

    sps[BUILD_CMD].add_argument(
        "-q", "--requirements", action="store_true",
        help="Show the build requirements for the specified asset and exit.")

    sps[BUILD_CMD].add_argument(
        "-r", "--recipe", required=False, default=None, type=str,
        help="Provide a recipe to use.")

    # add 'genome' argument to many commands
    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD, GETSEQ_CMD, TAG_CMD, ID_CMD]:
        # genome is not required for listing actions
        sps[cmd].add_argument(
            "-g", "--genome", required=cmd in GETSEQ_CMD,
            help="Reference assembly ID, e.g. mm10.")

    for cmd in LIST_REMOTE_CMD, LIST_LOCAL_CMD:
        sps[cmd].add_argument("-g", "--genome", required=False, type=str,
                              nargs="*", help="Reference assembly ID, e.g. mm10.")

    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD, TAG_CMD, ID_CMD]:
        sps[cmd].add_argument(
            "asset_registry_paths", metavar="asset-registry-paths", type=str, nargs='+',
            help="One or more registry path strings that identify assets  (e.g. hg38/fasta or hg38/fasta:tag"
                 + (" or hg38/fasta.fai:tag)." if cmd == GET_ASSET_CMD else ")."))

    for cmd in [PULL_CMD, REMOVE_CMD, INSERT_CMD]:
        sps[cmd].add_argument(
            "-f", "--force", action="store_true",
            help="Do not prompt before action, approve upfront.")

    sps[PULL_CMD].add_argument(
        "-u", "--no-untar", action="store_true",
        help="Do not extract tarballs.")

    sps[INSERT_CMD].add_argument(
        "-p", "--path", required=True,
        help="Relative local path to asset.")

    sps[GETSEQ_CMD].add_argument(
        "-l", "--locus", required=True,
        help="Coordinates of desired sequence; e.g. 'chr1:50000-50200'.")

    sps[GET_ASSET_CMD].add_argument(
        "-e", "--check-exists", required=False, action="store_true",
        help="Whether the returned asset path should be checked for existence "
             "on disk.")

    group = sps[TAG_CMD].add_mutually_exclusive_group(required=True)

    group.add_argument(
        "-t", "--tag", type=str,
        help="Tag to assign to an asset.")

    group.add_argument(
        "-d", "--default", action="store_true",
        help="Set the selected asset tag as the default one.")

    sps[SUBSCRIBE_CMD].add_argument(
        "-r", "--reset", action="store_true",
        help="Overwrite the current list of server URLs.")

    for cmd in [SUBSCRIBE_CMD, UNSUBSCRIBE_CMD]:
        sps[cmd].add_argument(
            "-s", "--genome-server", nargs='+', required=True,
            help="One or more URLs to {action} the {key} attribute in config file.".
                format(action="add to" if cmd == SUBSCRIBE_CMD else "remove from", key=CFG_SERVERS_KEY))

    return parser


def parse_registry_path(path):
    return prp(path, defaults=[
        ("protocol", None),
        ("genome", None),
        ("asset", None),
        ("seek_key", None),
        ("tag", None)])


def copy_or_download_file(input_string, outfolder):
    """
    Given an input file, which can be a local file or a URL, and output folder,
    this downloads or copies the file into the output folder.

    :param str input_string: Can be either a URL or a path to a local file
    :param str outfolder: Where to store the result.
    :return str, str: output/result file and command
    """
    result_file = os.path.join(outfolder, os.path.basename(input_string))
    parts = ["wget -O", result_file, input_string] \
        if is_url(input_string) else ["cp", input_string, result_file]
    return result_file, " ".join(parts)


def convert_file(input_fasta, output_file, conversions):
    """
    Given an input file, output file, and a list of conversions, gives the appropriate output file.

    :param str output_file: Path to local output file you want to create
    :param dict conversions: A dictionary of shell commands to convert files of a given type.
    """
    form = {"INPUT": input_fasta, "OUTPUT": output_file}
    _, ext = os.path.splitext(input_fasta)
    if ext in conversions:
        return conversions[ext].format(**form)


def default_config_file():
    """
    Path to default compute environment settings file.

    :return str: Path to default compute settings file
    """
    return os.path.join(os.path.dirname(__file__), "refgenie.yaml")


def get_asset_vars(genome, asset_key, tag, outfolder, specific_args=None, specific_params=None, **kwargs):
    """
    Gives a dict with variables used to populate an asset path.
    """
    asset_outfolder = os.path.join(outfolder, asset_key, tag)
    asset_vars = {"genome": genome,
                  "asset": asset_key,
                  "tag": tag,
                  "asset_outfolder": asset_outfolder}
    if specific_args:
        asset_vars.update(specific_args)
    if specific_params:
        asset_vars.update(specific_params)
    asset_vars.update(**kwargs)
    return asset_vars


def refgenie_add(rgc, asset_dict, path, force=False):
    """
    Add an external asset to the config.
    File existence is checked and asset files are transferred to the selected
    tag subdirectory

    :param refgenconf.RefGenConf rgc: genome configuration object
    :param dict asset_dict: a single parsed registry path
    :param str path: the path provided by the user. Must be relative to the
        specific genome directory
    :param bool force: whether the replacement of a possibly existing asset
        should be forced
    """
    # remove the first directory from the provided path if it is the genome name
    path = os.path.join(*path.split(os.sep)[1:]) \
        if path.split(os.sep)[0] == asset_dict["genome"] else path
    tag = asset_dict["tag"] \
          or rgc.get_default_tag(asset_dict["genome"], asset_dict["asset"])
    outfolder = \
        os.path.abspath(os.path.join(rgc[CFG_FOLDER_KEY], asset_dict["genome"]))
    abs_asset_path = os.path.join(outfolder, path)
    if asset_dict["seek_key"] is None:
        # if seek_key is not specified we're about to move a directory to
        # the tag subdir
        tag_path = os.path.join(abs_asset_path, tag)
        from shutil import copytree as cp
    else:
        # if seek_key is specified we're about to move just a single file to
        # he tag subdir
        tag_path = os.path.join(os.path.dirname(abs_asset_path), tag)
        if not os.path.exists(tag_path):
            os.makedirs(tag_path)
        from shutil import copy2 as cp
    if os.path.exists(abs_asset_path):
        if not os.path.exists(tag_path):
            cp(abs_asset_path, tag_path)
        else:
            if not force and not \
                    query_yes_no("Path '{}' exists. Do you want to overwrite?".
                                         format(tag_path)):
                return False
            else:
                _remove(tag_path)
                cp(abs_asset_path, tag_path)
    else:
        raise OSError("Absolute path '{}' does not exist. "
                      "The provided path must be relative to: {}".
                      format(abs_asset_path, rgc[CFG_FOLDER_KEY]))
    rgc.make_writable()
    gat_bundle = [asset_dict["genome"], asset_dict["asset"], tag]
    td = {CFG_ASSET_PATH_KEY:
              path if os.path.isdir(abs_asset_path) else os.path.dirname(path)}
    rgc.update_tags(*gat_bundle, data=td)
    # seek_key points to the entire dir if not specified
    seek_key_value = os.path.basename(abs_asset_path) \
        if asset_dict["seek_key"] is not None else "."
    sk = {asset_dict["seek_key"] or asset_dict["asset"]: seek_key_value}
    rgc.update_seek_keys(*gat_bundle, keys=sk)
    rgc.set_default_pointer(asset_dict["genome"], asset_dict["asset"], tag)
    # a separate update_tags call since we want to use the get_asset method
    # that requires a complete asset entry in rgc
    td = {CFG_ASSET_CHECKSUM_KEY: get_dir_digest(_seek(rgc, *gat_bundle))}
    rgc.update_tags(*gat_bundle, data=td)
    # Write the updated refgenie genome configuration
    rgc.write()
    rgc.make_readonly()
    return True


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
    genome_dir = os.path.join(rgc[CFG_FOLDER_KEY], genome)
    if is_writable(genome_dir):
        output_file = os.path.join(genome_dir, "{}_sequence_digests.tsv".format(genome))
        with open(output_file, "w") as contents_file:
            wr = csv.writer(contents_file, delimiter="\t")
            for key, val in content_checksums.items():
                wr.writerow([key, val])
        _LOGGER.debug("sequence digests saved to: {}".format(output_file))
    else:
        _LOGGER.warning("Could not save the genome sequence digests. '{}' is not writable".format(genome_dir))


def refgenie_build(gencfg, genome, asset_list, recipe_name, args):
    """
    Runs the refgenie build recipe.

    :param str gencfg: path to the genome configuration file
    :param argparse.Namespace args: parsed command-line options/arguments
    """
    rgc = RefGenConf(filepath=gencfg, writable=False)
    specified_args = _parse_user_build_input(args.files)
    specified_params = _parse_user_build_input(args.params)

    if not hasattr(args, "outfolder") or not args.outfolder:
        # Default to genome_folder
        _LOGGER.debug("No outfolder provided, using genome config.")
        args.outfolder = rgc[CFG_FOLDER_KEY]

    _LOGGER.debug("Default config file: {}".format(default_config_file()))

    if args.config_file and not os.path.isfile(args.config_file):
        _LOGGER.debug("Config file path isn't a file: {}".
                      format(args.config_file))
        args.config_file = default_config_file()

    def build_asset(genome, asset_key, tag, build_pkg, genome_outfolder, specific_args, specific_params, **kwargs):
        """
        Builds assets with pypiper and updates a genome config file.

        This function actually run the build commands in a given build package,
        and then update the refgenie config file.

        :param str genome: The assembly key; e.g. 'mm10'.
        :param str asset_key: The unique asset identifier; e.g. 'bowtie2_index'
        :param dict build_pkg: A dict (see examples) specifying lists
            of required input_assets, commands to run, and outputs to register as
            assets.
        """

        log_outfolder = os.path.abspath(os.path.join(genome_outfolder, asset_key, tag, BUILD_STATS_DIR))
        _LOGGER.info("Saving outputs to:\n- content: {}\n- logs: {}".format(genome_outfolder, log_outfolder))
        if args.docker:
            # Set up some docker stuff
            if args.volumes:
                # TODO: is volumes list defined here?
                volumes = volumes.append(genome_outfolder)
            else:
                volumes = genome_outfolder

        if not _writeable(genome_outfolder):
            _LOGGER.error("Insufficient permissions to write to output folder: {}".
                          format(genome_outfolder))
            return

        pm = pypiper.PipelineManager(name="refgenie", outfolder=log_outfolder, args=args)
        tk = pypiper.NGSTk(pm=pm)
        if args.docker:
            pm.get_container(build_pkg[CONT], volumes)
        _LOGGER.debug("Asset build package: " + str(build_pkg))
        gat = [genome, asset_key, tag]  # create a bundle list to simplify calls below
        # collect variables required to populate the command templates
        asset_vars = get_asset_vars(genome, asset_key, tag, genome_outfolder, specific_args, specific_params, **kwargs)
        # populate command templates
        # prior to populating, remove any seek_key parts from the keys, since these are not supported by format method
        command_list_populated = [x.format(**{k.split(".")[0]: v for k, v in asset_vars.items()})
                                  for x in build_pkg[CMD_LST]]
        # create output directory
        tk.make_dir(asset_vars["asset_outfolder"])

        target = os.path.join(log_outfolder, TEMPLATE_TARGET.format(genome, asset_key, tag))
        # add target command
        command_list_populated.append("touch {target}".format(target=target))
        _LOGGER.debug("Command populated: '{}'".format(" ".join(command_list_populated)))
        try:
            # run build command
            signal.signal(signal.SIGINT, _handle_sigint(gat))
            pm.run(command_list_populated, target, container=pm.container)
        except pypiper.exceptions.SubprocessError:
            _LOGGER.error("asset '{}' build failed".format(asset_key))
            return False
        else:
            # save build recipe to the JSON-formatted file
            recipe_file_name = TEMPLATE_RECIPE_JSON.format(asset_key, tag)
            with open(os.path.join(log_outfolder, recipe_file_name), 'w') as outfile:
                json.dump(build_pkg, outfile)
            # update and write refgenie genome configuration
            with rgc as r:
                r.update_assets(*gat[0:2], data={CFG_ASSET_DESC_KEY: build_pkg[DESC]})
                r.update_tags(*gat, data={CFG_ASSET_PATH_KEY: asset_key})
                r.update_seek_keys(*gat, keys={k: v.format(**asset_vars) for k, v in build_pkg[ASSETS].items()})
                # in order to conveniently get the path to digest we update the tags metadata in two steps
                digest = get_dir_digest(r.get_asset(genome, asset_key, tag, enclosing_dir=True), pm)
                r.update_tags(*gat, data={CFG_ASSET_CHECKSUM_KEY: digest})
                _LOGGER.info("Asset digest: {}".format(digest))
                r.set_default_pointer(*gat)
        pm.stop_pipeline()
        return True

    for a in asset_list:
        asset_key = a["asset"]
        asset_tag = a["tag"] or rgc.get_default_tag(genome, a["asset"], use_existing=False)
        recipe_name = recipe_name or asset_key

        if recipe_name in asset_build_packages.keys():
            asset_build_package = _check_recipe(asset_build_packages[recipe_name])
            # handle user-requested parents for the required assets
            input_assets = {}
            parent_assets = []
            specified_asset_keys, specified_assets = None, None
            if args.assets is not None:
                parsed_parents_input = _parse_user_build_input(args.assets)
                specified_asset_keys, specified_assets = \
                    list(parsed_parents_input.keys()), list(parsed_parents_input.values())
                _LOGGER.debug("Custom assets requested: {}".format(args.assets))
            if not specified_asset_keys and isinstance(args.assets, list):
                _LOGGER.warning("Specified parent assets format is invalid. Using defaults.")
            for req_asset in asset_build_package[REQ_ASSETS]:
                req_asset_data = parse_registry_path(req_asset[KEY])
                # for each req asset see if non-default parents were requested
                if specified_asset_keys is not None and req_asset_data["asset"] in specified_asset_keys:
                    parent_data = \
                        parse_registry_path(specified_assets[specified_asset_keys.index(req_asset_data["asset"])])
                    g, a, t, s = parent_data["genome"], \
                                 parent_data["asset"], \
                                 parent_data["tag"] or rgc.get_default_tag(genome, parent_data["asset"]), \
                                 parent_data["seek_key"]
                else:  # if no custom parents requested for the req asset, use default one
                    default = parse_registry_path(req_asset[DEFAULT])
                    g, a, t, s = genome, default["asset"], \
                                 rgc.get_default_tag(genome, default["asset"]), \
                                 req_asset_data["seek_key"]
                parent_assets.append("{}/{}:{}".format(g, a, t))
                input_assets[req_asset[KEY]] = _seek(rgc, g, a, t, s)
            _LOGGER.debug("Using parents: {}".format(", ".join(parent_assets)))
            _LOGGER.debug("Provided files: {}".format(specified_args))
            _LOGGER.debug("Provided parameters: {}".format(specified_params))
            for required_file in asset_build_package[REQ_FILES]:
                if specified_args is None or required_file[KEY] not in specified_args.keys():
                    raise ValueError("Path to the '{x}' input ({desc}) is required, but not provided. "
                                     "Specify it with: --files {x}=/path/to/{x}_file"
                                     .format(x=required_file[KEY], desc=required_file[DESC]))
            for required_param in asset_build_package[REQ_PARAMS]:
                if specified_params is None:
                    specified_params = {}
                if required_param[KEY] not in specified_params.keys():
                    if required_param[DEFAULT] is None:
                        raise ValueError("Value for the parameter '{x}' ({desc}) is required, but not provided. "
                                         "Specify it with: --params {x}=value"
                                         .format(x=required_param[KEY], desc=required_param[DESC]))
                    else:
                        specified_params.update({required_param[KEY]: required_param[DEFAULT]})
            genome_outfolder = os.path.join(args.outfolder, genome)
            _LOGGER.info("Building '{}/{}:{}' using '{}' recipe".format(genome, asset_key, asset_tag, recipe_name))
            if recipe_name == 'fasta' and genome in rgc.genomes_list() \
                    and 'fasta' in rgc.list_assets_by_genome(genome):
                _LOGGER.warning("'{g}' genome is already initialized with other fasta asset ({g}/{a}:{t}). "
                                "It will be re-initialized.".format(g=genome, a=asset_key, t=asset_tag))
            if not build_asset(genome, asset_key, asset_tag, asset_build_package, genome_outfolder,
                               specified_args, specified_params, **input_assets):
                log_path = os.path.abspath(os.path.join(genome_outfolder, asset_key, asset_tag,
                                                        BUILD_STATS_DIR, ORI_LOG_NAME))
                _LOGGER.info("'{}/{}:{}' was not added to the config, but directory has been left in place. "
                             "See the log file for details: {}".format(genome, asset_key, asset_tag, log_path))
                return
            # If the recipe was a fasta, we init the genome
            if recipe_name == 'fasta':
                _LOGGER.info("Computing initial genome digest...")
                collection_checksum, content_checksums = \
                    fasta_checksum(_seek(rgc, genome, asset_key, asset_tag, "fasta"))
                _LOGGER.info("Initializing genome...")
                refgenie_initg(rgc, genome, content_checksums)
            _LOGGER.info("Finished building '{}' asset".format(asset_key))
            with rgc as r:
                # update asset relationships
                r.update_relatives_assets(genome, asset_key, asset_tag, parent_assets)  # adds parents
                for i in parent_assets:
                    parsed_parent = parse_registry_path(i)
                    # adds child (currently built asset) to the parent
                    r.update_relatives_assets(parsed_parent["genome"], parsed_parent["asset"], parsed_parent["tag"],
                                                   ["{}/{}:{}".format(genome, asset_key, asset_tag)], True)
                if args.genome_description is not None:
                    _LOGGER.debug("adding genome ({}) description: '{}'".format(genome, args.genome_description))
                    r.update_genomes(genome, {CFG_GENOME_DESC_KEY: args.genome_description})
                if args.tag_description is not None:
                    _LOGGER.debug("adding tag ({}/{}:{}) description: '{}'".format(genome, asset_key, asset_tag,
                                                                                   args.tag_description))
                    r.update_tags(genome, asset_key, asset_tag, {CFG_TAG_DESC_KEY: args.tag_description})
                if recipe_name == "fasta":
                    # to save config lock time when building fasta assets
                    # (genome initialization takes some time for large genomes) we repeat the
                    # conditional here for writing the computed genome digest
                    r.update_genomes(genome, data={CFG_CHECKSUM_KEY: collection_checksum})
        else:
            _raise_missing_recipe_error(recipe_name)


def _exec_list(rgc, remote, genome):
    if remote:
        pfx = "Remote"
        assemblies, assets = rgc.get_remote_data_str(genome=genome)
        recipes = None  # Not implemented
    else:
        pfx = "Local"
        assemblies, assets = rgc.get_local_data_str(genome=genome)
        # also get recipes
        recipes = ", ".join(sorted(list(asset_build_packages.keys())))
    return pfx, assemblies, assets, recipes


def perm_check_x(file_to_check, message_tag):
    """
    Check X_OK permission on a path, providing according messaging and bool val.

    :param str file_to_check: path to query for permission
    :param str message_tag: context for error message if check fails
    :return bool: os.access(path, X_OK) for the given path
    :raise ValueError: if there's no filepath to check for permission
    """
    if not file_to_check:
        msg = "You must provide a path to {}".format(message_tag)
        _LOGGER.error(msg)
        raise ValueError(msg)
    if not os.access(file_to_check, os.X_OK):
        _LOGGER.error("Insufficient permissions to write to {}: "
                      "{}".format(message_tag, file_to_check))
        return False
    return True


def main():
    """ Primary workflow """
    parser = logmuse.add_logging_options(build_argparser())
    args, remaining_args = parser.parse_known_args()
    global _LOGGER
    _LOGGER = logmuse.logger_via_cli(args, make_root=True)
    _LOGGER.debug("refgenie {}".format(__version__))
    _LOGGER.debug("Args: {}".format(args))

    if not args.command:
        parser.print_help()
        _LOGGER.error("No command given")
        sys.exit(1)

    gencfg = refgenconf.select_genome_config(filename=args.genome_config, check_exist=not args.command == INIT_CMD,
                                             on_missing=lambda fp: fp, strict_env=True)
    if gencfg is None:
        raise MissingGenomeConfigError(args.genome_config)
    _LOGGER.debug("Determined genome config: {}".format(gencfg))

    # From user input we want to construct a list of asset dicts, where each
    # asset has a genome name, asset name, and tag

    if "asset_registry_paths" in args and args.asset_registry_paths:
        _LOGGER.debug("Found registry_path: {}".format(args.asset_registry_paths))
        asset_list = [parse_registry_path(x) for x in args.asset_registry_paths]

        for a in asset_list:
            # every asset must have a genome, either provided via registry path
            # or the args.genome arg.
            if not a["genome"]:
                if args.genome:
                    a["genome"] = args.genome
                else:
                    _LOGGER.error("Provided asset registry path ({}/{}:{}) is invalid. See help for usage reference.".
                                  format(a["genome"], a["asset"], a["tag"]))
                    sys.exit(1)
            else:
                if args.genome and args.genome != a["genome"]:
                    _LOGGER.warn("Two different genomes specified for asset '{}'.".format(a["asset"]))

    else:
        if args.command in GENOME_ONLY_REQUIRED and not args.genome:
            parser.error("You must provide either a genome or a registry path")
            sys.exit(1)
        if args.command in ASSET_REQUIRED:
            parser.error("You must provide an asset registry path")
            sys.exit(1)

    if args.command == INIT_CMD:
        _LOGGER.debug("Initializing refgenie genome configuration")
        rgc = RefGenConf(entries=OrderedDict({
            CFG_VERSION_KEY: REQ_CFG_VERSION,
            CFG_FOLDER_KEY: os.path.dirname(os.path.abspath(gencfg)),
            CFG_SERVERS_KEY: args.genome_server or [DEFAULT_SERVER],
            CFG_GENOMES_KEY: None
        }))
        rgc.initialize_config_file(os.path.abspath(gencfg))

    elif args.command == BUILD_CMD:
        if not all([x["genome"] == asset_list[0]["genome"] for x in asset_list]):
            _LOGGER.error("Build can only build assets for one genome")
            sys.exit(1)
        recipe_name = None
        if args.recipe:
            if len(asset_list) > 1:
                _LOGGER.error("Recipes cannot be specified for multi-asset builds")
                sys.exit(1)
            recipe_name = args.recipe
        if args.requirements:
            for a in asset_list:
                recipe = recipe_name or a["asset"]
                if recipe not in asset_build_packages.keys():
                    _raise_missing_recipe_error(recipe)
                _LOGGER.info("'{}' recipe requirements: ".format(recipe))
                _make_asset_build_reqs(recipe)
            sys.exit(0)
        refgenie_build(gencfg, asset_list[0]["genome"], asset_list, recipe_name, args)

    elif args.command == GET_ASSET_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        check = args.check_exists if args.check_exists else None
        for a in asset_list:
            _LOGGER.debug("getting asset: '{}/{}.{}:{}'".
                          format(a["genome"], a["asset"], a["seek_key"], a["tag"]))
            print(rgc.seek(a["genome"], a["asset"], a["tag"], a["seek_key"],
                           strict_exists=check))
        return

    elif args.command == INSERT_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        if len(asset_list) > 1:
            raise NotImplementedError("Can only add 1 asset at a time")
        else:
            refgenie_add(rgc, asset_list[0], args.path, args.force)

    elif args.command == PULL_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        force = None if not args.force else True
        outdir = rgc[CFG_FOLDER_KEY]
        if not os.path.exists(outdir):
            raise MissingFolderError(outdir)
        target = _key_to_name(CFG_FOLDER_KEY)
        if not perm_check_x(outdir, target):
            return
        if not _single_folder_writeable(outdir):
            _LOGGER.error("Insufficient permissions to write to {}: {}".
                          format(target, outdir))
            return

        for a in asset_list:
            rgc.pull(a["genome"], a["asset"], a["tag"],
                     unpack=not args.no_untar, force=force)

    elif args.command in [LIST_LOCAL_CMD, LIST_REMOTE_CMD]:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        if args.command == LIST_REMOTE_CMD:
            num_servers = 0
            # Keep all servers so that child updates maintain server list
            server_list = rgc[CFG_SERVERS_KEY]
            bad_servers = []
            for server_url in rgc[CFG_SERVERS_KEY]:
                num_servers += 1
                try:
                    rgc[CFG_SERVERS_KEY] = server_url
                    pfx, genomes, assets, recipes = _exec_list(rgc, args.command == LIST_REMOTE_CMD, args.genome)
                    if assets is None and genomes is None:
                        continue
                    _LOGGER.info("{} genomes: {}".format(pfx, genomes))
                    if args.command != LIST_REMOTE_CMD:  # Not implemented yet
                        _LOGGER.info("{} recipes: {}".format(pfx, recipes))
                    _LOGGER.info("{} assets:\n{}\n".format(pfx, assets))
                except (DownloadJsonError, ConnectionError):
                    bad_servers.append(server_url)
                    continue
            if num_servers >= len(server_list) and bad_servers:
                _LOGGER.error("Could not list assets from the following server(s): {}".format(bad_servers))
            # Restore original server list, even when we couldn't find assets on a server
            rgc[CFG_SERVERS_KEY] = server_list
        else:  # Only check local assets once
            _LOGGER.info("Server subscriptions: {}".format(", ".join(rgc[CFG_SERVERS_KEY])))
            pfx, genomes, assets, recipes = _exec_list(rgc, args.command == LIST_REMOTE_CMD, args.genome)
            _LOGGER.info("{} genomes: {}".format(pfx, genomes))
            if args.command != LIST_REMOTE_CMD:  # Not implemented yet
                _LOGGER.info("{} recipes: {}".format(pfx, recipes))
            _LOGGER.info("{} assets:\n{}".format(pfx, assets))

    elif args.command == GETSEQ_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        rgc.getseq(rgc, args.genome, args.locus)

    elif args.command == REMOVE_CMD:
        force = args.force
        rgc = RefGenConf(filepath=gencfg)
        for a in asset_list:
            a["tag"] = a["tag"] or rgc.get_default_tag(a["genome"], a["asset"],
                                                       use_existing=False)
            _LOGGER.debug("Determined tag for removal: {}".format(a["tag"]))
            if a["seek_key"] is not None:
                raise NotImplementedError("You can't remove a specific seek_key.")
            bundle = [a["genome"], a["asset"], a["tag"]]
            try:
                if not rgc.is_asset_complete(*bundle):
                    with rgc as r:
                        r.cfg_remove_assets(*bundle)
                    _LOGGER.info("Removed an incomplete asset '{}/{}:{}'".
                                 format(*bundle))
                    return
            except (KeyError, MissingAssetError, MissingGenomeError):
                _LOGGER.info("Asset '{}/{}:{}' does not exist".format(*bundle))
                return
        if len(asset_list) > 1:
            if not query_yes_no("Are you sure you want to remove {} assets?".
                                        format(len(asset_list))):
                _LOGGER.info("Action aborted by the user")
                return
            force = True
        for a in asset_list:
            rgc.remove(genome=a["genome"], asset=a["asset"], tag=a["tag"],
                       force=force)

    elif args.command == TAG_CMD:
        rgc = RefGenConf(filepath=gencfg)
        if len(asset_list) > 1:
            raise NotImplementedError("Can only tag 1 asset at a time")
        if args.default:
            # set the default tag and exit
            with rgc as r:
                r.set_default_pointer(a["genome"], a["asset"], a["tag"], True)
            sys.exit(0)
        rgc.tag(a["genome"], a["asset"], a["tag"], args.tag)

    elif args.command == ID_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        if len(asset_list) == 1:
            g, a = asset_list[0]["genome"], asset_list[0]["asset"]
            t = asset_list[0]["tag"] or rgc.get_default_tag(g, a)
            print(rgc.id(g, a, t))
            return
        for asset in asset_list:
            g, a = asset["genome"], asset["asset"]
            t = asset["tag"] or rgc.get_default_tag(g, a)
            print("{}/{}:{},".format(g, a, t) + rgc.id(g, a, t))
        return
    elif args.command == SUBSCRIBE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        rgc.subscribe(urls=args.genome_server, reset=args.reset)
        return
    elif args.command == UNSUBSCRIBE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        rgc.unsubscribe(urls=args.genome_server)
        return


def _entity_dir_removal_log(directory, entity_class, asset_dict, removed_entities):
    """
    Message and save removed entity data

    :param str directory: removed dir
    :param str entity_class: class of the entity
    :param dict asset_dict: selected genome/asset:tag combination
    :param list removed_entities: list of the removed entities to append to
    """
    subclass = "asset" if entity_class == "genome" else "tag"
    if os.path.basename(directory) == asset_dict[entity_class]:
        _LOGGER.info("Last {sub} for {ec} '{en}' has been removed, removing {ec} directory".
                     format(sub=subclass, ec=entity_class, en=asset_dict[entity_class]))
        removed_entities.append(_remove(directory))
    else:
        _LOGGER.debug("Didn't remove '{}' since it does not match the {} name: {}".
                      format(directory, entity_class, asset_dict[entity_class]))


def _remove(path):
    """
    remove asset if it is a dir or a file

    :param str path: path to the entity to remove, either a file or a dir
    :return str: removed path
    """
    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        rmtree(path)
    else:
        raise ValueError("path '{}' is neither a file nor a dir.".format(path))
    return path


def _key_to_name(k):
    return k.replace("_", " ")


def _single_folder_writeable(d):
    return os.access(d, os.W_OK) and os.access(d, os.X_OK)


def _writeable(outdir, strict_exists=False):
    outdir = outdir or "."
    if os.path.exists(outdir):
        return _single_folder_writeable(outdir)
    elif strict_exists:
        raise MissingFolderError(outdir)
    return _writeable(os.path.dirname(outdir), strict_exists)


def _make_asset_build_reqs(asset):
    """
    Prepare requirements and inputs lists and display it

    :params str asset: name of the asset
    """
    def _format_reqs(req_list):
        """

        :param list[dict] req_list:
        :return list[str]:
        """
        templ = "\t{} ({})"
        return [templ.format(req[KEY], req[DESC]) if DEFAULT not in req
                else (templ + "; default: {}").format(req[KEY], req[DESC], req[DEFAULT]) for req in req_list]

    reqs_list = []
    if asset_build_packages[asset][REQ_FILES]:
        reqs_list.append("- files:\n{}".format("\n".join(_format_reqs(asset_build_packages[asset][REQ_FILES]))))
    if asset_build_packages[asset][REQ_ASSETS]:
        reqs_list.append("- assets:\n{}".format("\n".join(_format_reqs(asset_build_packages[asset][REQ_ASSETS]))))
    if asset_build_packages[asset][REQ_PARAMS]:
        reqs_list.append("- params:\n{}".format("\n".join(_format_reqs(asset_build_packages[asset][REQ_PARAMS]))))
    _LOGGER.info("\n".join(reqs_list))


def get_dir_digest(path, pm=None):
    """
    Generate a MD5 digest that reflects just the contents of the files in the selected directory.

    :param str path: path to the directory to digest
    :param pypiper.PipelineManager pm: a pipeline object, optional. The subprocess module will be used if not provided
    :return str: a digest, e.g. a3c46f201a3ce7831d85cf4a125aa334
    """
    if not is_command_callable("md5sum"):
        raise OSError("md5sum command line tool is required for asset digest calculation. \n"
                      "Install and try again, e.g on macOS: 'brew install md5sha1sum'")
    cmd = "cd {}; find . -type f -not -path './" + BUILD_STATS_DIR + \
          "*' -exec md5sum {{}} \; | sort -k 2 | awk '{{print $1}}' | md5sum"
    if isinstance(pm, pypiper.PipelineManager):
        x = pm.checkprint(cmd.format(path))
    else:
        try:
            from subprocess import check_output
            x = check_output(cmd.format(path), shell=True).decode("utf-8")
        except Exception as e:
            _LOGGER.warning("{}: could not calculate digest for '{}'".format(e.__class__.__name__, path))
            return
    return str(sub(r'\W+', '', x))  # strips non-alphanumeric


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


def _parse_user_build_input(input):
    """
    Parse user input specification. Used in build for specific parents and input parsing.

    :param Iterable[Iterable[str], ...] input: user command line input,
        formatted as follows: [[fasta=txt, test=txt], ...]
    :return dict: mapping of keys, which are input names and values
    """
    lst = []
    for i in input or []:
        lst.extend(i)
    return {x.split("=")[0]: x.split("=")[1] for x in lst if "=" in x} if lst is not None else lst


def _raise_missing_recipe_error(recipe):
    """
    Raise an error for a missing recipe, when one is requested

    :param str recipe: recipe name
    :raise MissingRecipeError: always
    """
    raise MissingRecipeError("Recipe '{}' not found. Available recipes: {}".
                             format(recipe, ", ".join(list(asset_build_packages.keys()))))


def _check_recipe(recipe):
    """
    Check whether there are any key name clashes in the recipe requirements
    and raise an error if there are

    :param dict recipe: asset_build_package
    :raise ValueError: if any key names are duplicated
    """
    req_keys = []
    for req in [REQ_PARAMS, REQ_ASSETS, REQ_FILES]:
        req_keys.extend([req_dict[KEY] for req_dict in recipe[req]])
    unique = []
    for k in req_keys:
        if k not in unique:
            unique.append(k)
        else:
            raise ValueError("The recipe contains a duplicated requirement"
                             " key '{}', which is not permitted.".format(k))
    return recipe


def _seek(rgc, genome_name, asset_name, tag_name=None,
          seek_key=None, enclosing_dir=False):
    """
    Strict seek. Most use cases in this package require file existence
     check in seek. This function makes it easier
    """
    return rgc.seek(genome_name=genome_name,
                    asset_name=asset_name,
                    tag_name=tag_name,
                    seek_key=seek_key,
                    enclosing_dir=enclosing_dir,
                    strict_exists=True)
