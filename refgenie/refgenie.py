#!/usr/bin/env python

from argparse import SUPPRESS
from collections import OrderedDict
from shutil import rmtree
from re import sub
import os
import sys
import csv
import signal
import json

import pyfaidx

from ._version import __version__
from .exceptions import MissingGenomeConfigError, MissingFolderError
from .asset_build_packages import *

import logmuse
import pypiper
import refgenconf
from refgenconf import RefGenConf, MissingAssetError, MissingGenomeError, MissingRecipeError
from refgenconf.const import *
from ubiquerg import is_url, query_yes_no, parse_registry_path as prp, VersionInHelpParser, is_command_callable
from ubiquerg.system import is_writable
import yacman

# from refget import fasta_checksum
from .refget import fasta_checksum

_LOGGER = None

BUILD_CMD = "build"
INIT_CMD = "init"
PULL_CMD = "pull"
LIST_LOCAL_CMD = "list"
LIST_REMOTE_CMD = "listr"
GET_ASSET_CMD = "seek"
INSERT_CMD = "add"
REMOVE_CMD = "remove"
GETSEQ_CMD = "getseq"
TAG_CMD = "tag"

GENOME_ONLY_REQUIRED = [REMOVE_CMD, GETSEQ_CMD]

# For each asset we assume a genome is also required
ASSET_REQUIRED = [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, TAG_CMD]

BUILD_SPECIFIC_ARGS = ('fasta', 'ensembl_gtf', 'gencode_gtf', 'gff', 'context', 'refgene', 'dbnsfp')


def build_argparser():
    """
    Builds argument parser.

    :return argparse.ArgumentParser
    """

    banner = "%(prog)s - reference genome asset manager"
    additional_description = "\nhttps://refgenie.databio.org"

    parser = VersionInHelpParser(
        version=__version__,
        description=banner,
        epilog=additional_description)

    subparsers = parser.add_subparsers(dest="command")

    def add_subparser(cmd, description):
        return subparsers.add_parser(
            cmd, description=description, help=description)

    subparser_messages = {
        INIT_CMD: "Initialize a genome configuration.",
        LIST_LOCAL_CMD: "List available local assets.",
        LIST_REMOTE_CMD: "List available remote assets.",
        PULL_CMD: "Download assets.",
        BUILD_CMD: "Build genome assets.",
        GET_ASSET_CMD: "Get the path to a local asset.",
        INSERT_CMD: "Add local asset to the config file.",
        REMOVE_CMD: "Remove a local asset.",
        GETSEQ_CMD: "Get sequences from a genome.",
        TAG_CMD: "Tag an asset."
    }

    sps = {}
    for cmd, desc in subparser_messages.items():
        sps[cmd] = add_subparser(cmd, desc)
        # It's required for init
        sps[cmd].add_argument(
            '-c', '--genome-config', required=(cmd == INIT_CMD), dest="genome_config",
            help="Path to local genome configuration file. Optional if {} environment variable is set."
                .format(", ".join(refgenconf.CFG_ENV_VARS)))

    sps[INIT_CMD].add_argument('-s', '--genome-server', default=DEFAULT_SERVER,
                               help="URL to use for the genome_server attribute in config file. Default: {}"
                               .format(DEFAULT_SERVER))
    sps[BUILD_CMD] = pypiper.add_pypiper_args(
        sps[BUILD_CMD], groups=None, args=["recover", "config", "new-start"])

    # Add any arguments specific to subcommands.

    sps[BUILD_CMD].add_argument(
        '--tag-description', required=False, default=None, type=str,
        help="Add tag level description (e.g. built with version 0.3.2)")

    sps[BUILD_CMD].add_argument(
        '--genome-description', required=False, default=None, type=str,
        help="Add genome level description (e.g. The mouse mitochondrial genome, released in Dec 2013)")

    sps[BUILD_CMD].add_argument(
        "-d", "--docker", action="store_true", help="Run all commands in the refgenie docker container.")

    sps[BUILD_CMD].add_argument(
        '-t', '--tags', nargs="+", required=False, default=None,
        help='Override the default tags of the parent assets (e.g. asset:tag).')

    sps[BUILD_CMD].add_argument(
        '-v', '--volumes', nargs="+", required=False, default=None,
        help='If using docker, also mount these folders as volumes.')

    sps[BUILD_CMD].add_argument(
        '-o', '--outfolder', dest='outfolder', required=False, default=None,
        help='Override the default path to genomes folder, which is the '
             'genome_folder attribute in the genome configuration file.')

    sps[BUILD_CMD].add_argument(
        "-r", "--requirements", action="store_true",
        help="Show the build requirements for the specified asset and exit.")

    # add 'genome' argument to many commands
    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD,
                LIST_REMOTE_CMD, LIST_LOCAL_CMD, GETSEQ_CMD, TAG_CMD]:
        # genome is not required for listing actions
        sps[cmd].add_argument(
            "-g", "--genome", required=cmd in GETSEQ_CMD,
            help="Reference assembly ID, e.g. mm10")

    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD, TAG_CMD]:
        sps[cmd].add_argument(
            "asset_registry_paths", metavar="asset-registry-paths", type=str, nargs='+',
            help="One or more registry path strings that identify assets  (e.g. hg38/fasta or hg38/fasta:tag"
                 + (" or hg38/fasta.fai:tag)" if cmd == GET_ASSET_CMD else ")"))

    sps[PULL_CMD].add_argument(
        "-u", "--no-untar", action="store_true",
        help="Do not extract tarballs.")

    sps[INSERT_CMD].add_argument(
        "-p", "--path", required=True,
        help="Relative local path to asset.")

    sps[GETSEQ_CMD].add_argument(
        "-l", "--locus", required=True,
        help="Coordinates of desired sequence; e.g. 'chr1:50000-50200'.")

    group = sps[TAG_CMD].add_mutually_exclusive_group(required=True)

    group.add_argument(
        "-t", "--tag", type=str,
        help="Tag to assign to an asset.")

    group.add_argument(
        "-d", "--default", action="store_true",
        help="Set the selected asset tag as the default one.")

    # Finally, arguments to the build command to give the files needed to do
    # the building. These should eventually move to a more flexible system that
    # doesn't require them to be hard-coded here in order to be recognized

    for arg in BUILD_SPECIFIC_ARGS:
        sps[BUILD_CMD].add_argument(
            "--{arg}".format(arg=arg), required=False, help=SUPPRESS)

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


def get_asset_vars(genome, asset_key, tag, outfolder, specific_args=None, **kwargs):
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
    asset_vars.update(**kwargs)
    return asset_vars


def refgenie_add(rgc, asset_dict, path):
    """
    Add an external asset to the config.
    File existence is checked and asset files are transferred to the selected tag subdirectory

    :param refgenconf.RefGenConf rgc: genome configuration object
    :param dict asset_dict: a single parsed registry path
    :param str path: the path provided by the user. Must be relative to the specific genome directory
    """
    # remove the first directory from the provided path if it is the genome name
    path = os.path.join(*path.split(os.sep)[1:]) if path.split(os.sep)[0] == asset_dict["genome"] else path
    tag = asset_dict["tag"] or rgc.get_default_tag(asset_dict["genome"], asset_dict["asset"])
    outfolder = os.path.abspath(os.path.join(rgc.genome_folder, asset_dict["genome"]))
    abs_asset_path = os.path.join(outfolder, path)
    if asset_dict["seek_key"] is None:
        # if seek_key is not specified we're about to move a directory to the tag subdir
        tag_path = os.path.join(abs_asset_path, tag)
        from shutil import copytree as cp
    else:
        # if seek_key is not specified we're about to move just a single file to the tag subdir
        tag_path = os.path.join(os.path.dirname(abs_asset_path), tag)
        if not os.path.exists(tag_path):
            os.makedirs(tag_path)
        from shutil import copy2 as cp
    if os.path.exists(abs_asset_path):
        if not os.path.exists(tag_path):
            cp(abs_asset_path, tag_path)
        else:
            if not query_yes_no("Path '{}' exists. Do you want to overwrite?".format(tag_path)):
                return False
            else:
                _remove(tag_path)
                cp(abs_asset_path, tag_path)
    else:
        raise OSError("Absolute path '{}' does not exist. The provided path must be relative to: {}".
                      format(abs_asset_path, rgc.genome_folder))
    gat_bundle = [asset_dict["genome"], asset_dict["asset"], tag]
    rgc.update_tags(*gat_bundle,
                    data={CFG_ASSET_PATH_KEY: path if os.path.isdir(abs_asset_path) else os.path.dirname(path)})
    # seek_key points to the entire dir if not specified
    seek_key_value = os.path.basename(abs_asset_path) if asset_dict["seek_key"] is not None else "."
    rgc.update_seek_keys(*gat_bundle, keys={asset_dict["seek_key"] or asset_dict["asset"]: seek_key_value})
    rgc.set_default_pointer(asset_dict["genome"], asset_dict["asset"], tag)
    # a separate update_tags call since we want to use the get_asset method that requires a complete asset entry in rgc
    rgc.update_tags(*gat_bundle, data={CFG_ASSET_CHECKSUM_KEY: get_dir_digest(rgc.get_asset(*gat_bundle))})
    # Write the updated refgenie genome configuration
    rgc.write()
    return True


def refgenie_initg(rgc, genome, collection_checksum, content_checksums):
    """
    Initializing a genome means adding `collection_checksum` attributes in the
    genome config file. This should perhaps be a function in refgenconf, but not
    a CLI-hook. Also adds `content_checksums` tsv file (should be a recipe cmd?).

    This function updates the provided RefGenConf object with the
    genome(collection)-level checksum and saves the individual checksums to a
    TSV file in the fasta asset directory.

    :param refgenconf.RefGenConf rgc: genome configuration object
    :param str genome: name of the genome
    :param str collection_checksum: genome checksum
    :param dict content_checksums: checksums of individual content_checksums, e.g. chromosomes
    """
    rgc.update_genomes(genome, data={CFG_CHECKSUM_KEY: collection_checksum})
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


def refgenie_build(gencfg, genome, asset_list, args):
    """
    Runs the refgenie build recipe.

    :param str gencfg: path to the genome configuration file
    :param argparse.Namespace args: parsed command-line options/arguments
    """
    rgc = RefGenConf(filepath=gencfg, writable=False)
    specific_args = {k: getattr(args, k) for k in BUILD_SPECIFIC_ARGS}

    if not hasattr(args, "outfolder") or not args.outfolder:
        # Default to genome_folder
        _LOGGER.debug("No outfolder provided, using genome config.")
        args.outfolder = rgc.genome_folder

    _LOGGER.debug("Default config file: {}".format(default_config_file()))

    if args.config_file and not os.path.isfile(args.config_file):
        _LOGGER.debug("Config file path isn't a file: {}".
                      format(args.config_file))
        args.config_file = default_config_file()

    def build_asset(genome, asset_key, tag, build_pkg, genome_outfolder, specific_args, **kwargs):
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

        log_outfolder = os.path.abspath(os.path.join(genome_outfolder, asset_key, tag, "_refgenie_build"))
        _LOGGER.info("Output content: {}; logs: {}".format(genome_outfolder, log_outfolder))
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
        asset_vars = get_asset_vars(*gat, genome_outfolder, specific_args, **kwargs)
        # populate command templates
        command_list_populated = [x.format(**asset_vars) for x in build_pkg[CMD_LST]]
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
            rgc_rw = RefGenConf(filepath=gencfg, writable=True, wait_max=600)
            rgc_rw.update_assets(*gat[0:2], data={CFG_ASSET_DESC_KEY: build_pkg[DESC]})
            rgc_rw.update_tags(*gat, data={CFG_ASSET_PATH_KEY: asset_key})
            rgc_rw.update_seek_keys(*gat, keys={k: v.format(**asset_vars) for k, v in build_pkg[ASSETS].items()})
            # in order to conveniently get the path to digest we update the tags metadata in two steps
            digest = get_dir_digest(rgc_rw.get_asset(genome, asset_key, tag, enclosing_dir=True), pm)
            rgc_rw.update_tags(*gat, data={CFG_ASSET_CHECKSUM_KEY: digest})
            _LOGGER.info("Asset digest: {}".format(digest))
            rgc_rw.set_default_pointer(*gat)
            rgc_rw.write()
        pm.stop_pipeline()
        
        return True

    for a in asset_list:
        asset_key = a["asset"]
        asset_tag = a["tag"] or rgc.get_default_tag(genome, a["asset"], use_existing=False)

        if asset_key in asset_build_packages.keys():
            asset_build_package = asset_build_packages[asset_key]
            _LOGGER.debug(specific_args)
            # handle user-requested tags for the required assets
            input_assets = {}
            parent_assets = []
            parent_tags = args.tags
            selected_parent_tags = [p["item"] for p in [prp(x) for x in parent_tags]] if isinstance(parent_tags, list) \
                else None  # if tags specified construct asset_package.asset names
            for req_asset in asset_build_package[REQ_ASSETS]:
                req_asset_data = prp(req_asset)
                # for each req asset see if non-default tag was requested
                if selected_parent_tags is not None and req_asset_data["item"] in selected_parent_tags:
                    # if requested, add the path to the asset to the dictionary
                    parent_data = prp(parent_tags[asset_build_package[REQ_ASSETS].index(req_asset)])
                    if parent_data["tag"] is None:
                        raise ValueError("Parent asset tag was not specified. "
                                         "To use the default one, skip the -t/--tags option.")
                    input_assets[parent_data["item"]] = rgc.get_asset(genome, parent_data["item"], parent_data["tag"],
                                                                      req_asset_data["subitem"])
                    parent_assets.append("{}:{}".format(parent_data["item"], parent_data["tag"]))
                else:  # if no tag was requested for the req asset, use one tagged with default
                    default = prp(req_asset)
                    dafault_tag = rgc.get_default_tag(genome, default["item"])
                    input_assets[default["item"]] = rgc.get_asset(genome, default["item"], dafault_tag,
                                                                  default["subitem"])
                    parent_assets.append("{}:{}".format(default["item"], dafault_tag))
            _LOGGER.debug("parents: {}".format(", ".join(parent_assets)))
            _LOGGER.info("Inputs required to build '{}': {}".format(asset_key, ", ".join(asset_build_package[REQ_IN])))
            for required_input in asset_build_package[REQ_IN]:
                if not specific_args[required_input]:
                    raise ValueError(
                        "Argument '{}' is required to build asset '{}', but not provided".format(required_input,
                                                                                                 asset_key))
            error_req_template = "Asset '{}' is required to build asset '{}', but not provided"
            for ra in asset_build_package[REQ_ASSETS]:
                try:
                    req_data = prp(ra)
                    if not rgc.get_asset(genome, req_data["item"], req_data["tag"], req_data["subitem"]):
                        raise ValueError(error_req_template.format(ra, asset_key))
                except refgenconf.exceptions.MissingGenomeError:
                    raise ValueError(error_req_template.format(ra, asset_key))
            _LOGGER.info("Building asset '{}'".format(asset_key))
            genome_outfolder = os.path.join(args.outfolder, genome)
            if not build_asset(genome, asset_key, asset_tag, asset_build_package, genome_outfolder,
                               specific_args, **input_assets):
                log_path = os.path.abspath(os.path.join(genome_outfolder, asset_key, asset_tag,
                                                        "_refgenie_build", "refgenie_log.md"))
                _LOGGER.info("'{}/{}:{}' was not added to the config, but directory has been left in place. "
                             "See the log file for details: {}".format(genome, asset_key, asset_tag, log_path))
                return
            rgc_rw = RefGenConf(filepath=gencfg, writable=True, wait_max=600)  # create object for writing
            # If the asset was a fasta, we init the genome
            if asset_key == 'fasta':
                _LOGGER.info("Computing initial genome digest...")
                collection_checksum, content_checksums = \
                    fasta_checksum(rgc_rw.get_asset(genome, "fasta", asset_tag, "fasta"))
                _LOGGER.info("Initializing genome...")
                refgenie_initg(rgc_rw, genome, collection_checksum, content_checksums)
            _LOGGER.info("Finished building asset '{}'".format(asset_key))
            # update asset relationships
            rgc_rw.update_relatives_assets(genome, asset_key, asset_tag, parent_assets)  # adds parents
            for i in parent_assets:
                parsed_parent = prp(i)
                rgc_rw.update_relatives_assets(genome, parsed_parent["item"], parsed_parent["tag"],
                                            ["{}:{}".format(asset_key, asset_tag)], True)  # adds children
            if args.genome_description is not None:
                _LOGGER.debug("adding genome ({}) description: '{}'".format(genome, args.genome_description))
                rgc_rw.update_genomes(genome, {CFG_GENOME_DESC_KEY: args.genome_description})
            if args.tag_description is not None:
                _LOGGER.debug("adding tag ({}/{}:{}) description: '{}'".format(genome, asset_key, asset_tag,
                                                                               args.tag_description))
                rgc_rw.update_tags(genome, asset_key, asset_tag, {CFG_TAG_DESC_KEY: args.tag_description})
            rgc_rw.write()
            del rgc_rw
        else:
            raise MissingRecipeError("There is no '{}' recipe defined".format(asset_key))


def refgenie_init(genome_config_path, genome_server=DEFAULT_SERVER, config_version=REQ_CFG_VERSION):
    """
    Initialize a genome config file.

    :param str genome_config_path: path to genome configuration file to
        create/initialize
    :param str genome_server: URL for a server
    :param str config_version: config version name, e.g. 0.2
    """

    # Set up default
    rgc = RefGenConf(entries=OrderedDict({
        CFG_VERSION_KEY: config_version,
        CFG_FOLDER_KEY: os.path.dirname(os.path.abspath(genome_config_path)),
        CFG_SERVER_KEY: genome_server,
        CFG_GENOMES_KEY: None
    }))

    _LOGGER.debug("RGC: {}".format(rgc))

    if genome_config_path and not os.path.exists(genome_config_path):
        rgc.write(genome_config_path)
        _LOGGER.info("Wrote new refgenie genome configuration file: {}".format(genome_config_path))
    else:
        _LOGGER.warning("Can't initialize, file exists: {} ".format(genome_config_path))


def refgenie_getseq(rgc, genome, locus):
    """
    Something like the refget protocol.
    """

    fa = pyfaidx.Fasta(rgc.get_asset(genome, "fasta"))
    locus_split = locus.split(":")

    if len(locus_split) > 1:
        start, end = locus_split[1].split("-")
        _LOGGER.debug("chr: '{}', start: '{}', end: '{}'".format(locus_split[0], start, end))
        print(fa[locus_split[0]][int(start):int(end)])
    else:
        print(fa[locus_split[0]])


def _exec_list(rgc, remote, genome):
    if remote:
        pfx = "Remote"
        assemblies, assets = rgc.list_remote(genome=genome)
        recipes = None  # Not implemented
    else:
        pfx = "Local"
        assemblies, assets = rgc.list_local(genome=genome)
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

    gencfg = yacman.select_config(args.genome_config, CFG_ENV_VARS, check_exist=not args.command == INIT_CMD,
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
        _LOGGER.info("Initializing refgenie genome configuration")
        _writeable(os.path.dirname(gencfg), strict_exists=True)
        refgenie_init(gencfg, args.genome_server)
        sys.exit(0)

    if args.command == BUILD_CMD:
        if not all([x["genome"] == asset_list[0]["genome"] for x in asset_list]):
            _LOGGER.error("Build can only build assets from one genome")
            sys.exit(1)
        if args.requirements:
            if a["asset"] not in asset_build_packages.keys():
                _LOGGER.error("Recipe does not exist for asset '{}'".format(a["asset"]))
                sys.exit(1)
            _LOGGER.info("'{}/{}' build requirements: ".format(a["genome"], a["asset"]))
            _make_asset_build_reqs(a["asset"])
            sys.exit(0)
        refgenie_build(gencfg, asset_list[0]["genome"], asset_list, args)

    elif args.command == GET_ASSET_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        for a in asset_list:
            _LOGGER.debug("getting asset: '{}/{}.{}:{}'".format(a["genome"], a["asset"], a["seek_key"], a["tag"]))
            print(rgc.get_asset(a["genome"], a["asset"], a["tag"], a["seek_key"]))
        return

    elif args.command == INSERT_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=True)  # genome cfg will be updated, create object in RW mode mode
        if len(asset_list) > 1:
            raise NotImplementedError("Can only add 1 asset at a time")
        else:
            refgenie_add(rgc, asset_list[0], args.path)

    elif args.command == PULL_CMD:
        rgc = RefGenConf(filepath=gencfg)
        outdir = rgc[CFG_FOLDER_KEY]
        if not os.path.exists(outdir):
            raise MissingFolderError(outdir)
        target = _key_to_name(CFG_FOLDER_KEY)
        if not perm_check_x(outdir, target):
            return
        if not _single_folder_writeable(outdir):
            _LOGGER.error("Insufficient permissions to write to {}: {}".format(target, outdir))
            return

        for a in asset_list:
            gat, archive_data = rgc.pull_asset(a["genome"], a["asset"], a["tag"], unpack=not args.no_untar)
            if archive_data is not None:
                rgc_rw = RefGenConf(filepath=gencfg, writable=True)
                [rgc_rw.chk_digest_update_child(a["genome"], x, "{}:{}".format(gat[1], gat[2]))
                 for x in archive_data[CFG_ASSET_PARENTS_KEY] if CFG_ASSET_PARENTS_KEY in archive_data]
                rgc_rw.update_tags(*gat,
                                   data={attr: archive_data[attr] for attr in ATTRS_COPY_PULL if attr in archive_data})
                rgc_rw.set_default_pointer(*gat)
                rgc_rw.write()
                del rgc_rw

    elif args.command in [LIST_LOCAL_CMD, LIST_REMOTE_CMD]:
        rgc = RefGenConf(filepath=gencfg, writable=False)
        pfx, genomes, assets, recipes = _exec_list(rgc, args.command == LIST_REMOTE_CMD, args.genome)
        _LOGGER.info("{} genomes: {}".format(pfx, genomes))
        if args.command != LIST_REMOTE_CMD:  # Not implemented yet
            _LOGGER.info("{} recipes: {}".format(pfx, recipes))
        _LOGGER.info("{} assets:\n{}".format(pfx, assets))

    elif args.command == GETSEQ_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False)  # genome cfg will not be updated, create object in read-only mode
        refgenie_getseq(rgc, args.genome, args.locus)

    elif args.command == REMOVE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=True)  # genome cfg will be updated, create object in RW mode mode
        for a in asset_list:
            a["tag"] = a["tag"] or rgc.get_default_tag(a["genome"], a["asset"], use_existing=False)
            _LOGGER.debug("Determined tag for removal: {}".format(a["tag"]))
            if a["seek_key"] is not None:
                raise NotImplementedError("You can't remove a specific seek_key.")
            bundle = [a["genome"], a["asset"], a["tag"]]
            try:
                if not rgc.is_asset_complete(*bundle):
                    rgc.remove_assets(*bundle).write()
                    _LOGGER.info("Removed an incomplete asset '{}/{}:{}'".format(*bundle))
                    return
                else:
                    rgc.get_asset(*bundle, enclosing_dir=True)
            except (KeyError, MissingAssetError, MissingGenomeError):
                _LOGGER.info("Asset '{}/{}:{}' does not exist".format(*bundle))
                return
        if len(asset_list) > 1:
            if not query_yes_no("Are you sure you want to remove {} assets?".format(len(asset_list))):
                _LOGGER.info("Action aborted by the user")
                return
        else:
            a = asset_list[0]
            bundle = [a["genome"], a["asset"], a["tag"]]
            if not query_yes_no("Remove '{}/{}:{}'?".format(*bundle)):
                _LOGGER.info("Action aborted by the user")
                return
        removed = []
        for a in asset_list:
            bundle = [a["genome"], a["asset"], a["tag"]]
            asset_path = rgc.get_asset(*bundle, enclosing_dir=True)
            if os.path.exists(asset_path):
                removed.append(_remove(asset_path))
                rgc.remove_assets(*bundle).write()
            try:
                rgc[CFG_GENOMES_KEY][a["genome"]][CFG_ASSETS_KEY][a["asset"]]
            except (KeyError, TypeError):
                asset_dir = os.path.abspath(os.path.join(asset_path, os.path.pardir))
                _entity_dir_removal_log(asset_dir, "asset", a, removed)
                try:
                    rgc[CFG_GENOMES_KEY][a["genome"]][CFG_ASSETS_KEY]
                except (KeyError, TypeError):
                    genome_dir = os.path.abspath(os.path.join(asset_dir, os.path.pardir))
                    _entity_dir_removal_log(genome_dir, "genome", a, removed)
                    try:
                        del rgc[CFG_GENOMES_KEY][a["genome"]]
                        rgc.write()
                    except (KeyError, TypeError):
                        _LOGGER.debug("Could not remove genome '{}' from the config; it does not exist".
                                      format(a["genome"]))
            else:
                rgc.write()
        _LOGGER.info("Successfully removed entities:\n- {}".format("\n- ".join(removed)))

    elif args.command == TAG_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=True)  # genome cfg will be updated, create object in RW mode mode
        if len(asset_list) > 1:
            raise NotImplementedError("Can only tag 1 asset at a time")
        if args.default:
            # set the default tag and exit
            rgc.set_default_pointer(a["genome"], a["asset"], a["tag"], force=True)
            rgc.write()
            sys.exit(0)
        ori_path = rgc.get_asset(a["genome"], a["asset"], a["tag"], enclosing_dir=True)
        new_path = os.path.abspath(os.path.join(ori_path, os.pardir, args.tag))
        if not rgc.tag_asset(a["genome"], a["asset"], a["tag"], args.tag):  # tagging in the RefGenConf object
            sys.exit(0)
        try:
            if os.path.exists(new_path):
                _remove(new_path)
            os.rename(ori_path, new_path)  # tagging in the directory
        except FileNotFoundError:
            _LOGGER.warning("Could not rename original asset tag directory '{}' to the new one '{}'".
                            format(ori_path, new_path))
        else:
            rgc.remove_assets(a["genome"], a["asset"], a["tag"])
            _LOGGER.debug("Asset '{}/{}' tagged with '{}' has been removed from the genome config".
                          format(a["genome"], a["asset"], a["tag"]))
            _LOGGER.debug("Original asset has been moved from '{}' to '{}'".format(ori_path, new_path))
        rgc.write()


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

    def _format_req(req):
        """
        Format the asset requirements.

        Some of them specify seek_keys named as the assets, but we dont wanto to display that

        :param str req: requirement str
        :return str: formatted requirement
        """
        reqs = req.split(".")
        assert len(reqs) == 2, ValueError("Length of requirement '{}' after splitting is not 2. Specified requirement "
                                          "is invalid, it should be formatted as follows: 'asset.seek_key'".format(req))
        return reqs if reqs[1] != reqs[0] else reqs[1]

    reqs_list = []
    if asset_build_packages[asset][REQ_IN]:
        reqs_list.append("- arguments: {}".format(", ".join(asset_build_packages[asset][REQ_IN])))
    if asset_build_packages[asset][REQ_ASSETS]:
        reqs_list.append("- assets: {}".
                         format(", ".join(_format_req(r) for r in asset_build_packages[asset][REQ_ASSETS])))
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
    cmd = "cd {}; find . -type f -not -path './_refgenie_build*' " \
          "-exec md5sum {{}} \; |" \
          " sort -k 2 | awk '{{print $1}}' | " \
          "md5sum".format(path)
    if isinstance(pm, pypiper.PipelineManager):
        x = pm.checkprint(cmd)
    else:
        try:
            from subprocess import check_output
            x = check_output(cmd, shell=True).decode("utf-8")
        except Exception as e:
            _LOGGER.warning("{}: could not calculate digest for '{}'".format(e.__class__.__name__, path))
            return
    return sub(r'\W+', '', x)  # strips non-alphanumeric


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


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        _LOGGER.info("Program canceled by user!")
        sys.exit(1)
