#!/usr/bin/env python

from argparse import ArgumentParser, SUPPRESS
from collections import OrderedDict
from shutil import rmtree
import os
import re
import sys
import csv

import hashlib
import pyfaidx


from ._version import __version__
from .exceptions import MissingGenomeConfigError, MissingFolderError
from .asset_build_packages import *

import logmuse
import pypiper
import refgenconf
from refgenconf import RefGenConf, MissingAssetError, MissingGenomeError
from refgenconf.const import *
from ubiquerg import is_url, query_yes_no, parse_registry_path as prp, VersionInHelpParser
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
SETDEFAULT_CMD = "setdefault"

GENOME_ONLY_REQUIRED = [REMOVE_CMD, GETSEQ_CMD]

# For each asset we assume a genome is also required
ASSET_REQUIRED = [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, SETDEFAULT_CMD]

BUILD_SPECIFIC_ARGS = ('fasta', 'gtf', 'gff', 'context', 'refgene')

# This establishes the API with the server
refgenie_server_api = {
    "list_available_genomes": "/genomes",
    'list_assets_by_genome': "/genome/{genome}",
    'download_asset': "/asset/{genome}/{asset}",
}


def build_argparser():
    """
    Builds argument parser.

    :return argparse.ArgumentParser
    """

    banner = "%(prog)s - builds and manages reference genome assemblies"
    additional_description = "\nhttps://refgenie.databio.org"

    parser = VersionInHelpParser(
        description=banner,
        epilog=additional_description,
        version=__version__)

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
        GETSEQ_CMD: "Get sequences from a genome",
        SETDEFAULT_CMD: "Assign a 'default' tag to an asset"
    }

    sps = {}
    for cmd, desc in subparser_messages.items():
        sps[cmd] = add_subparser(cmd, desc)
        # It's required for init
        sps[cmd].add_argument(
            '-c', '--genome-config', required=(cmd == INIT_CMD), dest="genome_config",
            help="Path to local genome configuration file.")

    sps[INIT_CMD].add_argument('-s', '--genome-server', default=DEFAULT_SERVER,
                               help="URL to use for the genome_server attribute in config file."
                                    " Defaults : {}".format(DEFAULT_SERVER))
    sps[BUILD_CMD] = pypiper.add_pypiper_args(
        sps[BUILD_CMD], groups=None, args=["recover", "config", "new-start"])

    # Add any arguments specific to subcommands.

    sps[BUILD_CMD].add_argument(
        "-d", "--docker", action="store_true",
        help="Run all commands in the refgenie docker container.")

    sps[BUILD_CMD].add_argument(
        '-v', '--volumes', nargs="+", required=False, default=None,
        help='If using docker, also mount these folders as volumes')

    sps[BUILD_CMD].add_argument(
        '-o', '--outfolder', dest='outfolder', required=False, default=None,
        help='Override the default path to genomes folder, which is the '
             'genome_folder attribute in the genome configuration file.')

    # add 'genome' argument to many commands
    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD,
                 LIST_REMOTE_CMD, LIST_LOCAL_CMD, GETSEQ_CMD, SETDEFAULT_CMD]:
        # genome is not required for listing actions
        sps[cmd].add_argument(
            "-g", "--genome", required=cmd in (GETSEQ_CMD),
            help="Reference assembly ID, e.g. mm10")

    # add 'asset' argument to many commands
    # for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD]:
    #     sps[cmd].add_argument(
    #         "-a", "--asset", required=False, nargs='+',
    #         help="Name of one or more assets (keys in genome config file)")

    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, REMOVE_CMD, SETDEFAULT_CMD]:
        sps[cmd].add_argument(
            "asset_registry_paths", metavar="asset-registry-paths", type=str, nargs='+',
            help="One or more registry path strings that identify assets"
            " (e.g. hg38/bowtie2_index:1.0.0)")

    sps[PULL_CMD].add_argument(
        "-u", "--no-untar", action="store_true",
        help="Do not extract tarballs.")

    sps[INSERT_CMD].add_argument(
        "-p", "--path", required=True,
        help="Relative local path to asset")

    sps[GETSEQ_CMD].add_argument(
        "-l", "--locus", required=True,
        help="Coordinates to retrieve sequence for; such has 'chr1:50000-50200'.")


    # Finally, arguments to the build command to give the files needed to do
    # the building. These should eventually move to a more flexible system that
    # doesn't require them to be hard-coded here in order to be recognized

    for arg in BUILD_SPECIFIC_ARGS:
        sps[BUILD_CMD].add_argument(
            "--{arg}".format(arg=arg), required=False, help=SUPPRESS)

    return parser


def parse_registry_path(path):
    return prp(path, names=["protocol", "genome", "asset", "seek_key", "tag"])

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


def get_asset_vars(genome, asset_key, tag, outfolder, specific_args=None):
    """
    Gives a dict with variables used to populate an asset path.
    """
    asset_outfolder = os.path.join(outfolder, asset_key)
    asset_vars = {"genome": genome,
                  "asset": asset_key,
                  "tag": tag,
                  "asset_outfolder": asset_outfolder}
    if specific_args:
        asset_vars.update(specific_args)
    return asset_vars


def refgenie_add(rgc, asset_dict):
    outfolder = os.path.abspath(os.path.join(rgc.genome_folder, asset_dict["genome"]))
    rgc.update_assets(asset_dict["genome"],
                      asset_dict["asset"],
                      asset_dict["tag"],
                      {CFG_ASSET_PATH_KEY: args.path.format(**asset_dict)})

    # Write the updated refgenie genome configuration
    rgc.write()


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
    rgc.update_genomes(genome, {
            CFG_CHECKSUM_KEY: collection_checksum,
        })
    rgc.write()
    genome_dir = os.path.join(rgc[CFG_FOLDER_KEY], genome)
    if is_writable(genome_dir):
        output_file = os.path.join(genome_dir, "{}_sequence_checksums.tsv".format(genome))
        with open(output_file, "w") as contents_file:
            wr = csv.writer(contents_file, delimiter="\t")
            for key, val in content_checksums.items():
                wr.writerow([key, val])
        _LOGGER.debug("content checksums saved to: {}".format(output_file))
    else:
        _LOGGER.warning("Cound not save the genome content checksums hashes. "
                        "The directory '{}' os not writable".format(genome_dir))


def refgenie_build(rgc, genome, asset_list, args):
    """
    Runs the refgenie build recipe.

    :param refgenconf.RefGenConf rgc: genome configuration instance
    :param argparse.Namespace args: parsed command-line options/arguments
    """
    specific_args = {k: getattr(args, k) for k in BUILD_SPECIFIC_ARGS}

    if not hasattr(args, "outfolder") or not args.outfolder:
        # Default to genome_folder
        _LOGGER.debug("No outfolder provided, using genome config.")
        args.outfolder = rgc.genome_folder

    outfolder = os.path.abspath(os.path.join(args.outfolder, genome))
    if not _writeable(outfolder):
        _LOGGER.error("Insufficient permissions to write to output folder: {}".
                      format(outfolder))
        return

    _LOGGER.info("Output to: {} {} {}".format(genome, args.outfolder, outfolder))
    _LOGGER.debug("Default config file: {}".format(default_config_file()))

    if args.config_file and not os.path.isfile(args.config_file):
        _LOGGER.debug("Config file path isn't a file: {}".
                      format(args.config_file))
        args.config_file = default_config_file()

    def build_asset(genome, asset_key, tag, build_pkg, outfolder, specific_args):
        """
        Builds assets with pypiper and updates a genome config file.

        This function actually run the build commands in a given build package,
        and then update the refgenie config file.

        :param str genome: The assembly key; e.g. 'mm10'.
        :param str asset_key: The unique asset identifier; e.g. 'bowtie2_index'
        :param dict build_pkg: A dict (see examples) specifying lists
            of required inputs, commands to run, and outputs to register as
            assets.
        """
        _LOGGER.debug("Asset build package: " + str(build_pkg))
        asset_vars = get_asset_vars(genome, asset_key, tag, outfolder, specific_args)
        asset_outfolder = os.path.join(outfolder, asset_key)

        _LOGGER.debug(str([x.format(**asset_vars) for x in build_pkg[CMD_LST]]))

        tk.make_dir(asset_outfolder)
        target = os.path.join(asset_outfolder, "build_complete.flag")
        command_list_populated = [x.format(**asset_vars) for x in build_pkg[CMD_LST]]

        touch_target = "touch {target}".format(target=target)
        command_list_populated.append(touch_target)

        _LOGGER.debug("Command list populated: " + str(command_list_populated))

        pm.run(command_list_populated, target, container=pm.container)

        for asset in build_pkg[ASSETS].keys():
            rgc.update_assets(genome, asset, tag, {
                CFG_ASSET_PATH_KEY: build_pkg[ASSETS][asset][PTH].format(**asset_vars),
                CFG_ASSET_DESC_KEY: build_pkg[ASSETS][asset][ASSET_DESC]
            })

        # Write the updated refgenie genome configuration
        rgc.write()

    pm = pypiper.PipelineManager(name="refgenie", outfolder=outfolder, args=args)
    tk = pypiper.NGSTk(pm=pm)

    if args.docker:
        # Set up some docker stuff
        if args.volumes:
            volumes = volumes.append(outfolder)
        else:
            volumes = outfolder

    for a in asset_list:
        asset_key = a["asset"]
        asset_tag = a["tag"]

        if asset_key in asset_build_packages.keys():
            asset_build_package = asset_build_packages[asset_key]
            _LOGGER.debug(specific_args)
            required_inputs = ", ".join(asset_build_package[REQ_IN])
            _LOGGER.info("Inputs required to build '{}': {}".format(asset_key, required_inputs))
            for required_input in asset_build_package[REQ_IN]:
                if not specific_args[required_input]:
                    raise ValueError(
                        "Argument '{}' is required to build asset '{}', but not provided".format(required_input,
                                                                                                 asset_key))
            for required_asset in asset_build_package[REQ_ASSETS]:
                try:
                    if not rgc.get_asset(args.genome, required_asset):
                        raise ValueError(
                            "Asset '{}' is required to build asset '{}', but not provided".format(required_asset,
                                                                                                  asset_key))
                except refgenconf.exceptions.MissingGenomeError:
                    raise ValueError(
                        "Asset '{}' is required to build asset '{}', but not provided".format(required_asset,
                                                                                              asset_key))
            if args.docker:
                pm.get_container(asset_build_package[CONT], volumes)

            # If the asset is a fasta, we first init the genome
            if asset_key == 'fasta':
                _LOGGER.info("Initializing genome...")
                collection_checksum, content_checksums = fasta_checksum(specific_args["fasta"])
                if genome in rgc.genomes and CFG_CHECKSUM_KEY in rgc.genomes[genome]\
                        and collection_checksum != rgc.genomes[genome][CFG_CHECKSUM_KEY]:
                    _LOGGER.info("Checksum doesn't match")
                    return False
                refgenie_initg(rgc, genome, collection_checksum, content_checksums)
            build_asset(genome, asset_key, asset_tag, asset_build_package, outfolder, specific_args)
            _LOGGER.info("Finished building asset '{}'".format(asset_key))
        else:
            _LOGGER.warn("Recipe does not exist for asset '{}'".format(asset_key))

    pm.stop_pipeline()


def refgenie_init(genome_config_path, genome_server=DEFAULT_SERVER, config_version=REQ_CFG_VERSION):
    """
    Initialize a genome config file.

    :param str genome_config_path: path to genome configuration file to
        create/initialize
    :param str genome_server: URL for a server
    :param str config_version: config version name, e.g. 0.2
    """

    # Set up default
    rgc = RefGenConf(OrderedDict({
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
        recipes = ", ".join(list(asset_build_packages.keys()))

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
    _LOGGER = logmuse.logger_via_cli(args)
    logmuse.logger_via_cli(args, name=refgenconf.__name__)
    _LOGGER.debug("refgenie {}".format(__version__))
    _LOGGER.debug("Args: {}".format(args))

    if not args.command:
        parser.print_help()
        _LOGGER.error("No command given")
        sys.exit(1)

    gencfg = yacman.select_config(
        args.genome_config, CFG_ENV_VARS,
        check_exist=not args.command == INIT_CMD, on_missing=lambda fp: fp)
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
                    _LOGGER.error("Asset lacks a genome: {}".format(a["asset"]))
                    sys.exit(1)
            else:
                if args.genome and args.genome != a["genome"]:
                    _LOGGER.warn("Genome specified twice for asset '{}'.".format(
                        a["name"]))

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

    rgc = RefGenConf(gencfg)

    if args.command == BUILD_CMD:
        if not all([x["genome"] == asset_list[0]["genome"] for x in asset_list]):
            _LOGGER.error("Build can only build assets from one genome")
            sys.exit(1)
        refgenie_build(rgc, asset_list[0]["genome"], asset_list, args)

    elif args.command == GET_ASSET_CMD:
        for a in asset_list:
            _LOGGER.debug("getting asset: '{}/{}.{}:{}'".format(a["genome"], a["asset"], a["seek_key"], a["tag"]))
            print(rgc.get_asset(a["genome"], a["asset"], a["seek_key"], a["tag"]))
        return

    elif args.command == INSERT_CMD:
        if len(asset_list) > 1:
            raise NotImplementedError("Can only add 1 asset at a time")
        else:
            refgenie_add(rgc, asset_list[0])

    elif args.command == PULL_CMD:
        outdir = rgc[CFG_FOLDER_KEY]
        if not os.path.exists(outdir):
            raise MissingFolderError(outdir)
        target = _key_to_name(CFG_FOLDER_KEY)
        if not perm_check_x(outdir, target):
            return
        if not _single_folder_writeable(outdir):
            _LOGGER.error("Insufficient permissions to write to {}: "
                          "{}".format(target, outdir))
            return

        for a in asset_list:
            rgc.pull_asset(a["genome"], a["asset"], a["tag"], gencfg,
                           unpack=not args.no_untar)

    elif args.command in [LIST_LOCAL_CMD, LIST_REMOTE_CMD]:
        pfx, genomes, assets, recipes = _exec_list(rgc, args.command == LIST_REMOTE_CMD, args.genome)
        _LOGGER.info("{} genomes: {}".format(pfx, genomes))
        if args.command != LIST_REMOTE_CMD:  # Not implemented yet
            _LOGGER.info("{} recipes: {}".format(pfx, recipes))
        _LOGGER.info("{} assets:\n{}".format(pfx, assets))

    elif args.command == GETSEQ_CMD:
        refgenie_getseq(rgc, args.genome, args.locus)

    elif args.command == REMOVE_CMD:
        if len(asset_list) < 1:
            # No assets provided, must be all for this genome
            asset_list = rgc.list_assets_by_genome(args.genome)

        for asset in assets:
            try:
                rgc.get_asset(args.genome, asset)
            except (MissingAssetError, MissingGenomeError):
                _LOGGER.info("Asset {}/{} does not exist".format(args.genome, asset))
                return
        if len(assets) > 1:
            if not query_yes_no("Remove {} assets for genome: {}?".format(len(assets), args.genome)):
                return
        else:
            query_yes_no("Remove {}/{}?".format(args.genome, assets[0]))
        removed = []
        for asset in assets:
            asset_dir = os.path.abspath(os.path.join(rgc.get_asset(args.genome, asset), os.pardir))
            asset_archive = asset_dir + ".tar"
            if os.path.isfile(asset_archive):
                os.remove(asset_archive)
                removed.append(asset_archive)
            if os.path.isdir(asset_dir):
                rmtree(asset_dir)
                rgc.remove_assets(args.genome, asset).write()
                removed.append(asset_dir)
        _LOGGER.info("Successfully removed entities:\n- {}".format("\n- ".join(removed)))

    elif args.command == SETDEFAULT_CMD:
        for a in asset_list:
            rgc.set_default_asset(a["genome"], a["asset"], a["tag"])
            rgc.write()

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


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        _LOGGER.info("Program canceled by user!")
        sys.exit(1)
