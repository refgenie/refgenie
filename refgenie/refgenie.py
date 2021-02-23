#!/usr/bin/env python

from collections import OrderedDict
from shutil import rmtree
from re import sub
from requests import ConnectionError
from rich.console import Console

import os
import sys
import csv
import signal
import json

from ._version import __version__
from .exceptions import MissingGenomeConfigError, MissingFolderError
from .asset_build_packages import *
from .const import *

import logmuse
import pypiper
import refgenconf
from refgenconf import RefGenConf, MissingAssetError, MissingGenomeError, \
    MissingRecipeError, DownloadJsonError, get_dir_digest, upgrade_config, \
    __version__ as rgc_version, select_genome_config
from ubiquerg import is_url, query_yes_no, parse_registry_path as prp, \
    VersionInHelpParser, is_command_callable
from ubiquerg.system import is_writable
from yacman import UndefinedAliasError
from argparse import HelpFormatter

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
        version=f"{__version__} | refgenconf {rgc_version}",
        description=banner,
        epilog=additional_description)

    subparsers = parser.add_subparsers(dest="command")

    def add_subparser(cmd, msg, subparsers):
        return subparsers.add_parser(
            cmd, description=msg, help=msg,
            formatter_class=lambda prog: HelpFormatter(
                prog, max_help_position=40, width=90
            )
        )

    sps = {}
    for cmd, desc in SUBPARSER_MESSAGES.items():
        sps[cmd] = add_subparser(cmd, desc, subparsers)
        # alias is nested and alias subcommands require config path
        if cmd == ALIAS_CMD:
            continue
        # It's required for init
        sps[cmd].add_argument(
            '-c', '--genome-config', required=(cmd == INIT_CMD), dest="genome_config", metavar="C",
            help="Path to local genome configuration file. Optional if {} environment variable is set."
                .format(", ".join(CFG_ENV_VARS)))
        sps[cmd].add_argument(
            '--skip-read-lock', required=False, action="store_true",
            help="Whether the config file should not be locked for reading")

    # upgrade: upgrade config and alter file structure to the target version
    sps[UPGRADE_CMD].add_argument('-v', '--target-version', required=True, metavar="V",
                                  help="Target config version for the upgrade.")
    sps[UPGRADE_CMD].add_argument('-f', '--force', action="store_true",
                                  help="Do not prompt before action, approve upfront.")

    sps[INIT_CMD].add_argument('-s', '--genome-server', nargs='+', default=DEFAULT_SERVER,
                               help="URL(s) to use for the {} attribute in config file. Default: {}."
                               .format(CFG_SERVERS_KEY, DEFAULT_SERVER))
    sps[INIT_CMD].add_argument('-f', '--genome-folder',
                               help="Absolute path to parent folder refgenie-managed assets.")
    sps[INIT_CMD].add_argument('-a', '--genome-archive-folder',
                               help="Absolute path to parent archive folder refgenie-managed assets; used by refgenieserver.")
    sps[INIT_CMD].add_argument('-b', '--genome-archive-config',
                               help="Absolute path to desired archive config file; used by refgenieserver.")
    sps[INIT_CMD].add_argument('-u', '--remote-url-base',
                               help="URL to use as an alternative, remote archive location; used by refgenieserver.")
    sps[INIT_CMD].add_argument('-j', '--settings-json',
                               help="Absolute path to a JSON file with the key "
                                    "value pairs to inialize the configuration "
                                    "file with. Overwritten by itemized specifications.")
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

    alias_subparser = sps[ALIAS_CMD]
    alias_subsubparsers = alias_subparser.add_subparsers(dest="subcommand")

    alias_sps = {}
    for cmd, desc in ALIAS_SUBPARSER_MESSAGES.items():
        alias_sps[cmd] = add_subparser(cmd, desc, alias_subsubparsers)
        alias_sps[cmd].add_argument(
            '-c', '--genome-config', required=False, dest="genome_config", metavar="C",
            help="Path to local genome configuration file. Optional if {} environment variable is set."
            .format(", ".join(CFG_ENV_VARS)))
        alias_sps[cmd].add_argument(
            '--skip-read-lock', required=False, action="store_true",
            help="Whether the config file should not be locked for reading")

    alias_sps[ALIAS_SET_CMD].add_argument(
        "-a", "--aliases", metavar="A", required=False, default=None, type=str,
        nargs="+", help="Aliases to set; single if the digest is to be retrieved from the server.")
    alias_sps[ALIAS_SET_CMD].add_argument(
        "-d", "--digest", metavar="D", required=False, type=str,
        help="Digest to set; leave out if the digest is to be retrieved from the server.")
    alias_sps[ALIAS_SET_CMD].add_argument(
        "-r", "--reset", action="store_true",
        help="Whether all the aliases should be removed prior to setting new ones.")
    alias_sps[ALIAS_SET_CMD].add_argument(
        "-f", "--force", action="store_true",
        help="Whether the action should be forced, if genome does not exist.")

    alias_sps[ALIAS_REMOVE_CMD].add_argument(
        "-a", "--aliases", metavar="A", required=False, default=None, type=str,
        nargs="+", help="Aliases to remove.")
    alias_sps[ALIAS_REMOVE_CMD].add_argument(
        "-d", "--digest", metavar="D", required=True, type=str,
        help="Digest to remove.")

    alias_sps[ALIAS_GET_CMD].add_argument(
        "-a", "--aliases", metavar="A", required=False, type=str, nargs="+",
        help="Aliases to get the digests for.")

    sps[COMPARE_CMD].add_argument("genome1", metavar="GENOME1", type=str, nargs=1,
                                  help="First genome for compatibility check.")
    sps[COMPARE_CMD].add_argument("genome2", metavar="GENOME2", type=str, nargs=1,
                                  help="Second genome for compatibility check.")
    sps[COMPARE_CMD].add_argument("-e", "--no-explanation", action="store_true",
                                  help="Do not print compatibility code explanation.")

    # add 'genome' argument to many commands
    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD, GETSEQ_CMD, TAG_CMD, ID_CMD]:
        # genome is not required for listing actions
        sps[cmd].add_argument(
            "-g", "--genome", required=cmd in GETSEQ_CMD,  metavar="G",
            help="Reference assembly ID, e.g. mm10.")

    for cmd in LIST_REMOTE_CMD, LIST_LOCAL_CMD:
        sps[cmd].add_argument("-g", "--genome", required=False, type=str, metavar="G",
                              nargs="*", help="Reference assembly ID, e.g. mm10.")

    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD, TAG_CMD, ID_CMD]:
        sps[cmd].add_argument(
            "asset_registry_paths", metavar="asset-registry-paths", type=str, nargs='+',
            help="One or more registry path strings that identify assets  (e.g. hg38/fasta or hg38/fasta:tag"
                 + (" or hg38/fasta.fai:tag)." if cmd == GET_ASSET_CMD else ")."))

    sps[LIST_LOCAL_CMD].add_argument("-r", "--recipes", action="store_true",
                                     help="List available recipes.")

    for cmd in [REMOVE_CMD, INSERT_CMD]:
        sps[cmd].add_argument(
            "-f", "--force", action="store_true",
            help="Do not prompt before action, approve upfront.")

    sps[REMOVE_CMD].add_argument(
        "-a", "--aliases", action="store_true",
        help="Remove the genome alias if last asset for that genome is removed.")
    force_group = sps[PULL_CMD].add_argument_group(
        title="Prompt handling",
        description="These flags configure the pull prompt responses.")

    overwrite_group = force_group.add_mutually_exclusive_group()

    overwrite_group.add_argument("--no-overwrite", action="store_true",
                                 help="Do not overwrite if asset exists.")

    overwrite_group.add_argument("--force-overwrite", action="store_true",
                                 help="Overwrite if asset exists.")

    large_group = force_group.add_mutually_exclusive_group()

    large_group.add_argument("--no-large", action="store_true",
                             help="Do not pull archives over 5GB.")

    large_group.add_argument("--pull-large", action="store_true",
                             help="Pull any archive, regardless of its size.")

    force_group.add_argument("--size-cutoff", type=float, default=10, metavar="S",
                             help="Maximum archive file size to download with no confirmation required (in GB, default: 10)")

    force_group.add_argument("-b", "--batch", action="store_true",
                             help="Use batch mode: pull large archives, do no overwrite")

    sps[INSERT_CMD].add_argument(
        "-p", "--path", required=True, metavar="P",
        help="Relative local path to asset.")

    sps[INSERT_CMD].add_argument(
        "-s", "--seek-keys", required=False, type=str, metavar="S",
        help="""
        String representation of a JSON object with seek_keys, 
        e.g. '{"seek_key1": "file.txt"}'
        """)

    sps[GETSEQ_CMD].add_argument(
        "-l", "--locus", required=True,
        help="Coordinates of desired sequence; e.g. 'chr1:50000-50200'.")

    sps[GET_ASSET_CMD].add_argument(
        "-e", "--check-exists", required=False, action="store_true",
        help="Whether the returned asset path should be checked for existence on disk.")

    sps[TAG_CMD].add_argument(
        '-f', '--force', action="store_true",
        help="Do not prompt before action, approve upfront.")

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
        output_file = os.path.join(
            genome_dir, "{}_sequence_digests.tsv".format(genome))
        with open(output_file, "w") as contents_file:
            wr = csv.writer(contents_file, delimiter="\t")
            for key, val in content_checksums.items():
                wr.writerow([key, val])
        _LOGGER.debug("sequence digests saved to: {}".format(output_file))
    else:
        _LOGGER.warning(
            "Could not save the genome sequence digests. '{}' is not writable".format(genome_dir))


def refgenie_build(gencfg, genome, asset_list, recipe_name, args):
    """
    Runs the refgenie build recipe.

    :param str gencfg: path to the genome configuration file
    :param argparse.Namespace args: parsed command-line options/arguments
    """
    rgc = RefGenConf(filepath=gencfg, writable=False,
                     skip_read_lock=_skip_lock(args.skip_read_lock, gencfg))
    specified_args = _parse_user_build_input(args.files)
    specified_params = _parse_user_build_input(args.params)

    def _read_json_file(filepath):
        """
        Read a JSON file

        :param str filepath: path to the file to read
        :return dict: read data
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
        return data

    if recipe_name and os.path.isfile(recipe_name) and recipe_name.endswith(".json"):
        recipe_name = _read_json_file(filepath=recipe_name)

    if not hasattr(args, "outfolder") or not args.outfolder:
        # Default to genome_folder
        _LOGGER.debug("No outfolder provided, using genome config.")
        args.outfolder = rgc.data_dir

    _LOGGER.debug("Default config file: {}".format(default_config_file()))

    if args.config_file and not os.path.isfile(args.config_file):
        _LOGGER.debug("Config file path isn't a file: {}".
                      format(args.config_file))
        args.config_file = default_config_file()

    def _build_asset(genome, asset_key, tag, build_pkg, genome_outfolder, specific_args, specific_params, alias, **kwargs):
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

        log_outfolder = os.path.abspath(os.path.join(
            genome_outfolder, asset_key, tag, BUILD_STATS_DIR))
        _LOGGER.info(
            "Saving outputs to:\n- content: {}\n- logs: {}".format(genome_outfolder, log_outfolder))
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

        pm = pypiper.PipelineManager(
            name="refgenie", outfolder=log_outfolder, args=args)
        tk = pypiper.NGSTk(pm=pm)
        if args.docker:
            pm.get_container(build_pkg[CONT], volumes)
        _LOGGER.debug("Asset build package: " + str(build_pkg))
        # create a bundle list to simplify calls below
        gat = [genome, asset_key, tag]
        # collect variables required to populate the command templates
        asset_vars = get_asset_vars(
            genome, asset_key, tag, genome_outfolder, specific_args, specific_params, **kwargs)
        # populate command templates
        # prior to populating, remove any seek_key parts from the keys, since these are not supported by format method
        command_list_populated = [x.format(**{k.split(".")[0]: v for k, v in asset_vars.items()})
                                  for x in build_pkg[CMD_LST]]
        # create output directory
        tk.make_dir(asset_vars["asset_outfolder"])

        target = os.path.join(
            log_outfolder, TEMPLATE_TARGET.format(genome, asset_key, tag))
        # add target command
        command_list_populated.append("touch {target}".format(target=target))
        _LOGGER.debug("Command populated: '{}'".format(
            " ".join(command_list_populated)))
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
            # since the assets are always built to a standard dir structure, we
            # can just stitch a path together for asset digest calculation
            asset_dir = os.path.join(rgc.data_dir, *gat)
            if not os.path.exists(asset_dir):
                raise OSError("Could not compute asset digest. Path does not "
                              "exist: {}".format(asset_dir))
            digest = get_dir_digest(asset_dir)
            _LOGGER.info("Asset digest: {}".format(digest))
            # add updates to config file
            with rgc as r:
                if asset_key == "fasta":
                    r.update_genomes(genome, data={CFG_ALIASES_KEY: [alias]},
                                     force_digest=genome)
                r.update_assets(
                    *gat[0:2], data={CFG_ASSET_DESC_KEY: build_pkg[DESC]},
                    force_digest=genome)
                r.update_tags(
                    *gat, force_digest=genome,
                    data={CFG_ASSET_PATH_KEY: asset_key, CFG_ASSET_CHECKSUM_KEY: digest})
                r.update_seek_keys(
                    *gat, force_digest=genome,
                    keys={k: v.format(**asset_vars) for k, v in build_pkg[ASSETS].items()})
                r.set_default_pointer(*gat, force_digest=genome)
        pm.stop_pipeline()
        return True

    for a in asset_list:
        asset_key = a["asset"]
        asset_tag = a["tag"] or \
                    rgc.get_default_tag(genome, a["asset"], use_existing=False)
        recipe_name = recipe_name or asset_key

        if isinstance(recipe_name, dict) or \
                (isinstance(recipe_name, str)
                 and recipe_name in asset_build_packages.keys()):
            if isinstance(recipe_name, dict):
                _LOGGER.info("Using custom recipe: \n{}".format(recipe_name))
                asset_build_package = _check_recipe(recipe_name)
                recipe_name = asset_build_package["name"]
            else:
                asset_build_package = \
                    _check_recipe(asset_build_packages[recipe_name])
            # handle user-requested parents for the required assets
            input_assets = {}
            parent_assets = []
            specified_asset_keys, specified_assets = None, None
            if args.assets is not None:
                parsed_parents_input = _parse_user_build_input(args.assets)
                specified_asset_keys = list(parsed_parents_input.keys())
                specified_assets = list(parsed_parents_input.values())
                _LOGGER.debug(f"Custom assets requested: {args.assets}")
            if not specified_asset_keys and isinstance(args.assets, list):
                _LOGGER.warning(
                    "Specified parent assets format is invalid. Using defaults.")
            for req_asset in asset_build_package[REQ_ASSETS]:
                req_asset_data = parse_registry_path(req_asset[KEY])
                # for each req asset see if non-default parents were requested
                if specified_asset_keys is not None and req_asset_data["asset"] in specified_asset_keys:
                    parent_data = parse_registry_path(
                        specified_assets[specified_asset_keys.index(req_asset_data["asset"])])
                    g, a, t, s = parent_data["genome"], \
                        parent_data["asset"], \
                        parent_data["tag"] or rgc.get_default_tag(genome, parent_data["asset"]), \
                        parent_data["seek_key"]
                else:  # if no custom parents requested for the req asset, use default one
                    default = parse_registry_path(req_asset[DEFAULT])
                    g, a, t, s = genome, default["asset"], \
                        rgc.get_default_tag(genome, default["asset"]), \
                        req_asset_data["seek_key"]
                parent_assets.append(
                    "{}/{}:{}".format(rgc.get_genome_alias_digest(g, fallback=True), a, t))
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
            _LOGGER.info("Building '{}/{}:{}' using '{}' recipe".format(
                genome, asset_key, asset_tag, recipe_name))
            ori_genome = genome
            if recipe_name == 'fasta':
                if genome in rgc.genomes_list() and 'fasta' in rgc.list_assets_by_genome(genome):
                    pretag = rgc.get_default_tag(genome, "fasta")
                    _LOGGER.warning("'{g}' genome is already initialized with other fasta asset ({g}/{a}:{t})".
                                    format(g=genome, a=asset_key, t=pretag))
                    genome = rgc.get_genome_alias_digest(alias=genome, fallback=True)
                else:
                    # if the recipe is "fasta" we first initialiaze the genome, based on the provided path to the input FASTA file
                    genome, _ = rgc.initialize_genome(
                        fasta_path=specified_args["fasta"], alias=ori_genome, skip_alias_write=True)
            else:
                try:
                    genome = rgc.get_genome_alias_digest(genome, fallback=True)
                except UndefinedAliasError:
                    _LOGGER.error("Genome '{}' has not been initialized yet; "
                                  "no key found for this alias".format(genome))
                    return
            recipe_name = None
            genome_outfolder = os.path.join(args.outfolder, genome)
            if not _build_asset(genome, asset_key, asset_tag, asset_build_package, genome_outfolder,
                                specified_args, specified_params, ori_genome, **input_assets):
                log_path = os.path.abspath(os.path.join(genome_outfolder, asset_key, asset_tag,
                                                        BUILD_STATS_DIR, ORI_LOG_NAME))
                _LOGGER.info("'{}/{}:{}' was not added to the config, but directory has been left in place. "
                             "See the log file for details: {}".format(genome, asset_key, asset_tag, log_path))
                return
            _LOGGER.info("Finished building '{}' asset".format(asset_key))
            with rgc as r:
                # update asset relationships
                r.update_relatives_assets(
                    genome, asset_key, asset_tag, parent_assets)  # adds parents
                for i in parent_assets:
                    parsed_parent = parse_registry_path(i)
                    # adds child (currently built asset) to the parent
                    r.update_relatives_assets(parsed_parent["genome"], parsed_parent["asset"],
                                              parsed_parent["tag"], ["{}/{}:{}".format(genome, asset_key, asset_tag)], True)
                if args.genome_description is not None:
                    _LOGGER.debug("adding genome ({}) description: '{}'".format(
                        genome, args.genome_description))
                    r.update_genomes(
                        genome, {CFG_GENOME_DESC_KEY: args.genome_description})
                if args.tag_description is not None:
                    _LOGGER.debug("adding tag ({}/{}:{}) description: '{}'".
                                  format(genome, asset_key, asset_tag, args.tag_description))
                    r.update_tags(genome, asset_key, asset_tag, {
                                  CFG_TAG_DESC_KEY: args.tag_description})
            rgc._symlink_alias(genome, asset_key, asset_tag)
        else:
            _raise_missing_recipe_error(recipe_name)


def _exec_list(rgc, remote, genome):
    if remote:
        pfx = "Remote"
        # we use this func looping through the server urls and assigning a
        # single instance as the server for the object. That's why we can
        # access the data with [0] below
        assemblies, assets = \
            list(rgc.listr(genome=genome, as_str=True).values())[0]
        recipes = None  # Not implemented
    else:
        pfx = "Local"
        assemblies, assets = rgc.get_local_data_str(genome=genome)
        # also get recipes
        recipes = ", ".join(sorted(list(asset_build_packages.keys())))
    return pfx, assemblies, assets, recipes


def perm_check_x(file_to_check, message_tag="genome directory"):
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
        _LOGGER.error(
            "Insufficient permissions to write to {}: ".format(file_to_check))
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

    if args.command == ALIAS_CMD and not args.subcommand:
        parser.print_help()
        _LOGGER.error("No alias subcommand command given")
        sys.exit(1)

    gencfg = select_genome_config(
        filename=args.genome_config, check_exist=not args.command == INIT_CMD,
        on_missing=lambda fp: fp, strict_env=True)
    if gencfg is None:
        raise MissingGenomeConfigError(args.genome_config)
    _LOGGER.debug("Determined genome config: {}".format(gencfg))

    skip_read_lock = _skip_lock(args.skip_read_lock, gencfg)

    # From user input we want to construct a list of asset dicts, where each
    # asset has a genome name, asset name, and tag
    if "asset_registry_paths" in args and args.asset_registry_paths:
        _LOGGER.debug("Found registry_path: {}".format(
            args.asset_registry_paths))
        asset_list = [parse_registry_path(x)
                      for x in args.asset_registry_paths]

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
                    _LOGGER.warn(
                        "Two different genomes specified for asset '{}'.".format(a["asset"]))

    else:
        if args.command in GENOME_ONLY_REQUIRED and not args.genome:
            parser.error("You must provide either a genome or a registry path")
            sys.exit(1)
        if args.command in ASSET_REQUIRED:
            parser.error("You must provide an asset registry path")
            sys.exit(1)

    if args.command == INIT_CMD:
        _LOGGER.debug("Initializing refgenie genome configuration")
        entries = OrderedDict({
            CFG_VERSION_KEY: REQ_CFG_VERSION,
            CFG_FOLDER_KEY: os.path.dirname(os.path.abspath(gencfg)),
            CFG_SERVERS_KEY: args.genome_server or [DEFAULT_SERVER],
            CFG_GENOMES_KEY: None})
        if args.settings_json:
            if os.path.isfile(args.settings_json):
                with open(args.settings_json, 'r') as json_file:
                    data = json.load(json_file)
                entries.update(data)
            else:
                raise FileNotFoundError(
                    "JSON file with config init settings does not exist: {}".
                    format(args.settings_json))
        if args.genome_folder:
            entries.update({CFG_FOLDER_KEY: args.genome_folder})
        if args.remote_url_base:
            entries.update({CFG_REMOTE_URL_BASE_KEY: args.remote_url_base})
        if args.genome_archive_folder:
            entries.update({CFG_ARCHIVE_KEY: args.genome_archive_folder})
        if args.genome_archive_config:
            entries.update(
                {CFG_ARCHIVE_CONFIG_KEY: args.genome_archive_config})
        _LOGGER.debug("initializing with entries: {}".format(entries))
        rgc = RefGenConf(entries=entries, skip_read_lock=skip_read_lock)
        rgc.initialize_config_file(os.path.abspath(gencfg))

    elif args.command == BUILD_CMD:
        if not all([x["genome"] == asset_list[0]["genome"] for x in asset_list]):
            _LOGGER.error("Build can only build assets for one genome")
            sys.exit(1)
        recipe_name = None
        if args.recipe:
            if len(asset_list) > 1:
                _LOGGER.error(
                    "Recipes cannot be specified for multi-asset builds")
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
        refgenie_build(
            gencfg, asset_list[0]["genome"], asset_list, recipe_name, args)

    elif args.command == GET_ASSET_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
        check = args.check_exists if args.check_exists else None
        for a in asset_list:
            _LOGGER.debug("getting asset: '{}/{}.{}:{}'".
                          format(a["genome"], a["asset"], a["seek_key"], a["tag"]))
            print(rgc.seek(a["genome"], a["asset"], a["tag"], a["seek_key"],
                           strict_exists=check))
        return

    elif args.command == INSERT_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)

        if len(asset_list) > 1:
            raise NotImplementedError("Can only add 1 asset at a time")
        else:
            sk = args.seek_keys
            if sk:
                sk = json.loads(args.seek_keys)
            rgc.add(path=args.path, genome=asset_list[0]["genome"],
                    asset=asset_list[0]["asset"], tag=asset_list[0]["tag"],
                    seek_keys=sk, force=args.force)

    elif args.command == PULL_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)

        # existing assets overwriting
        if args.no_overwrite:
            force = False
        elif args.force_overwrite:
            force = True
        else:
            force = None
        # large archive pulling
        if args.no_large:
            force_large = False
        elif args.pull_large:
            force_large = True
        else:
            force_large = None
        # batch mode takes precedence over other choices
        if args.batch:
            force_large = True
            force = False

        outdir = rgc.data_dir
        if not os.path.exists(outdir):
            raise MissingFolderError(outdir)
        if not perm_check_x(outdir):
            return
        if not _single_folder_writeable(outdir):
            _LOGGER.error(
                "Insufficient permissions to write to: {}".format(outdir))
            return

        for a in asset_list:
            rgc.pull(a["genome"], a["asset"], a["tag"], force=force,
                     force_large=force_large, size_cutoff=args.size_cutoff)

    elif args.command in [LIST_LOCAL_CMD, LIST_REMOTE_CMD]:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
        console = Console()
        if args.command == LIST_REMOTE_CMD:
            num_servers = 0
            bad_servers = []
            for server_url in rgc[CFG_SERVERS_KEY]:
                num_servers += 1
                try:
                    table = rgc.get_asset_table(
                        genomes=args.genome, server_url=server_url)
                except (DownloadJsonError, ConnectionError):
                    bad_servers.append(server_url)
                    continue
                else:
                    console.print(table)
            if num_servers >= len(rgc[CFG_SERVERS_KEY]) and bad_servers:
                _LOGGER.error(
                    "Could not list assets from the following servers: {}".
                    format(bad_servers)
                )
        else:
            if args.recipes:
                print(", ".join(sorted(list(asset_build_packages.keys()))))
            else:
                console.print(rgc.get_asset_table(genomes=args.genome))

    elif args.command == GETSEQ_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
        print(rgc.getseq(args.genome, args.locus))

    elif args.command == REMOVE_CMD:
        force = args.force
        rgc = RefGenConf(filepath=gencfg,
                         skip_read_lock=skip_read_lock)
        for a in asset_list:
            a["tag"] = a["tag"] or rgc.get_default_tag(a["genome"], a["asset"],
                                                       use_existing=False)
            _LOGGER.debug("Determined tag for removal: {}".format(a["tag"]))
            if a["seek_key"] is not None:
                raise NotImplementedError(
                    "You can't remove a specific seek_key.")
            gat = {"genome": a["genome"], "asset": a["asset"], "tag": a["tag"]}
            try:
                if not rgc.is_asset_complete(**gat):
                    with rgc as r:
                        r.cfg_remove_assets(**gat, aliases=args.aliases)
                    _LOGGER.info("Removed an incomplete asset "
                                 "'{genome}/{asset}:{tag}'".format(*gat))
                    return
            except (KeyError, MissingAssetError, MissingGenomeError):
                _LOGGER.info("Asset '{genome}/{asset}:{tag}' does not exist"
                             .format(**gat))
                return
        if len(asset_list) > 1:
            if not query_yes_no("Are you sure you want to remove {} assets?".
                                format(len(asset_list))):
                _LOGGER.info("Action aborted by the user")
                return
            force = True
        for a in asset_list:
            rgc.remove(genome=a["genome"], asset=a["asset"],
                       tag=a["tag"], force=force)

    elif args.command == TAG_CMD:
        rgc = RefGenConf(filepath=gencfg,
                         skip_read_lock=skip_read_lock)
        if len(asset_list) > 1:
            raise NotImplementedError("Can only tag 1 asset at a time")
        if args.default:
            # set the default tag and exit
            with rgc as r:
                r.set_default_pointer(a["genome"], a["asset"], a["tag"], True)
            sys.exit(0)
        rgc.tag(a["genome"], a["asset"], a["tag"], args.tag, force=args.force)

    elif args.command == ID_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
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
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
        rgc.subscribe(urls=args.genome_server, reset=args.reset)
        return
    elif args.command == UNSUBSCRIBE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
        rgc.unsubscribe(urls=args.genome_server)
        return
    elif args.command == ALIAS_CMD:
        rgc = RefGenConf(filepath=gencfg, skip_read_lock=skip_read_lock)
        if args.subcommand == ALIAS_GET_CMD:
            if args.aliases is not None:
                for a in args.aliases:
                    print(rgc.get_genome_alias_digest(alias=a))
                return
            console = Console()
            console.print(rgc.genome_aliases_table)

        if args.subcommand == ALIAS_SET_CMD:
            rgc.set_genome_alias(digest=args.digest, genome=args.aliases,
                                 reset_digest=args.reset, create_genome=args.force)
            return
        elif args.subcommand == ALIAS_REMOVE_CMD:
            rgc.remove_genome_aliases(digest=args.digest, aliases=args.aliases)
            return

    elif args.command == COMPARE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False,
                         skip_read_lock=skip_read_lock)
        res = rgc.compare(args.genome1[0], args.genome2[0],
                          explain=not args.no_explanation)
        if args.no_explanation:
            print(res)

    elif args.command == UPGRADE_CMD:
        upgrade_config(target_version=args.target_version, filepath=gencfg,
                       force=args.force)


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
        reqs_list.append(
            "- files:\n{}".format("\n".join(_format_reqs(asset_build_packages[asset][REQ_FILES]))))
    if asset_build_packages[asset][REQ_ASSETS]:
        reqs_list.append(
            "- assets:\n{}".format("\n".join(_format_reqs(asset_build_packages[asset][REQ_ASSETS]))))
    if asset_build_packages[asset][REQ_PARAMS]:
        reqs_list.append(
            "- params:\n{}".format("\n".join(_format_reqs(asset_build_packages[asset][REQ_PARAMS]))))
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
            _LOGGER.warning("{}: could not calculate digest for '{}'".format(
                e.__class__.__name__, path))
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
    # experimental feature; recipe jsonschema validation
    from jsonschema import validate
    from yacman import load_yaml
    SCHEMA_SRC = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "schemas", "recipe_schema.yaml")
    if os.path.exists(SCHEMA_SRC):
        validate(recipe, load_yaml(filepath=SCHEMA_SRC))
        _LOGGER.info(
            "Recipe validated successfully against a schema: {}".format(SCHEMA_SRC))
    else:
        _LOGGER.warning("Recipe schema not found: {}".format(SCHEMA_SRC))
    # end of validation
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
    return rgc.seek_src(genome_name=genome_name,
                        asset_name=asset_name,
                        tag_name=tag_name,
                        seek_key=seek_key,
                        enclosing_dir=enclosing_dir,
                        strict_exists=True)


def _skip_lock(skip_arg, cfg):
    """
    If config read lock skip was not forced, check if dir is writable and set
    the default to the result

    :param bool skip_arg: argument selected on the CLI
    :param str cfg: path to the confjg
    :return bool: decision -- whether to skip the file lock for read
    """
    return is_writable(os.path.dirname(cfg)) if not skip_arg else True
