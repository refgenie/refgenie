#!/usr/bin/env python

from argparse import ArgumentParser, SUPPRESS
from collections import OrderedDict
import os
import re
import sys

from ._version import __version__
from .exceptions import MissingGenomeConfigError, MissingFolderError

import logmuse
import pypiper
import refgenconf
from refgenconf import RefGenConf
from refgenconf.const import *
from ubiquerg import is_url
import yacman

_LOGGER = None

BUILD_CMD = "build"
INIT_CMD = "init"
PULL_CMD = "pull"
LIST_LOCAL_CMD = "list"
LIST_REMOTE_CMD = "listr"
GET_ASSET_CMD = "seek"


BUILD_SPECIFIC_ARGS = ('fasta', 'gtf', 'context')

# This establishes the API with the server
refgenie_server_api = {
    "list_available_genomes": "/genomes",
    'list_assets_by_genome': "/genome/{genome}",
    'download_asset': "/asset/{genome}/{asset}",
}


class _VersionInHelpParser(ArgumentParser):
    def format_help(self):
        """ Add version information to help text. """
        return "version: {}\n".format(__version__) + \
               super(_VersionInHelpParser, self).format_help()


def build_argparser():
    """
    Builds argument parser.

    :return argparse.ArgumentParser
    """

    banner = "%(prog)s - builds and manages reference genome assemblies"
    additional_description = "\nhttps://refgenie.databio.org"

    parser = _VersionInHelpParser(
        description=banner,
        epilog=additional_description)

    parser.add_argument(
        "-V", "--version",
        action="version",
        version="%(prog)s {v}".format(v=__version__))

    subparsers = parser.add_subparsers(dest="command")

    def add_subparser(cmd, description):
        return subparsers.add_parser(
            cmd, description=description, help=description)

    subparser_messages = {
        INIT_CMD: "Initialize a genome configuration.",
        LIST_LOCAL_CMD: "List available local genomes.",
        LIST_REMOTE_CMD: "List available genomes and assets on server.",
        PULL_CMD: "Download assets.",
        BUILD_CMD: "Build genome assets",
        GET_ASSET_CMD: "Get the path to a local asset"
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
        '-o', '--outfolder', dest='outfolder', required=False, default=None,
        help='Override the default path to genomes folder, which is the '
             'genome_folder attribute in the genome configuration file.')

    for cmd in [PULL_CMD, GET_ASSET_CMD, BUILD_CMD]:
        sps[cmd].add_argument(
            "-g", "--genome", required=True,
            help="Reference assembly ID, e.g. mm10")
        sps[cmd].add_argument(
            "-a", "--asset", required=True, nargs='+',
            help="Name of one or more assets (keys in genome config file)")

    sps[PULL_CMD].add_argument(
        "-u", "--no-untar", action="store_true",
        help="Do not extract tarballs.")

    # Finally, arguments to the build command to give the files needed to do
    # the building. These should eventually move to a more flexible system that
    # doesn't require them to be hard-coded here in order to be recognized

    for arg in BUILD_SPECIFIC_ARGS:
        sps[BUILD_CMD].add_argument(
            "--{arg}".format(arg=arg), required=False, help=SUPPRESS)

    # sps[BUILD_CMD].add_argument(
    #     '--fasta', required=False, help=SUPPRESS)
    # help='Local path or URL to genome sequence file in .fa, .fa.gz, '
    #          'or .2bit format.'
    # sps[BUILD_CMD].add_argument(
    #        '--gtf', required=False, help=SUPPRESS)
    # help='Path to GTF gene annotation file.'


    return parser


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


def refgenie_build(rgc, args):
    """
    Runs the refgenie build recipe.
    
    :param refgenconf.RefGenConf rgc: genome configuration instance
    :param argparse.Namespace args: parsed command-line options/arguments
    """
     
    # Build specific args

    specific_args = {k: getattr(args,k) for k in BUILD_SPECIFIC_ARGS}


    if args.genome:
        genome = args.genome
    else:
        # This can probably be eliminated now that with flexible building
        genome = os.path.basename(args.input)
        # eliminate extensions to get canonical genome name.
        for strike in [".fasta.gz$", ".fa.gz$", ".fasta$", ".fa$", ".gz$", ".2bit$"]:
            genome = re.sub(strike, "", genome)

    _LOGGER.info("Using genome name: {}".format(genome))

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

    def path_data(root, c):
        return {"path": os.path.relpath(root, c.genome_folder)}



# This dict provides 'asset packages', which specify recipes (commands) to
# build assets. Each package can produce one or more assets, which are encoded
# as relative paths. The package name is often the same as the asset name but
# it need not be. 
    asset_build_packages = {
        "fasta": {
            "assets": {
                "fasta": "fasta/{genome}.fa",
                "fai": "fasta/{genome}.fa.fai",
                "chrom_sizes": "fasta/{genome}.chrom.sizes",
            },
            "required_inputs": ["fasta"],
            "command_list": [
                "cp {fasta} {asset_outfolder}/{genome}.fa",
                "samtools faidx {asset_outfolder}/{genome}.fa",
                "cut -f 1,2 {asset_outfolder}/{genome}.fa.fai > {asset_outfolder}/{genome}.chrom.sizes"
            ]
        },
        "bowtie2_index": {
            "assets": {
                "bowtie2_index": "bowtie2_index",
            },       
            "required_inputs": ["fasta"],
            "command_list": [
                "bowtie2-build {fasta} {asset_outfolder}/{genome}"
                ] 
        },
        "hisat2_index": {
            "assets": {
                "hisat2_index": "hisat2_index",
            },       
            "required_inputs": ["fasta"],
            "command_list": [
                "hisat2-build {fasta} {asset_outfolder}/{genome}"
                ] 
        },
        "bismark_bt2_index": {
            "assets": {
                "bismark_bt2_index": "bismark_bt2_index",
            },       
            "required_inputs": ["fasta"],
            "command_list": [
                "ln -sf ../{genome}.fa.gz {asset_outfolder}"
                "bismark_genome_preparation --bowtie2 {fasta} {asset_outfolder}/{genome}"
                ] 
        },
        "bismark_bt1_index": {
            "assets": {
                "bismark_bt1_index": "bismark_bt1_index",
            },       
            "required_inputs": ["fasta"],
            "command_list": [
                "ln -sf ../{genome}.fa.gz {asset_outfolder}"
                "bismark_genome_preparation {fasta} {asset_outfolder}/{genome}"
                ] 
        },  
        "kallisto_index": {
            "required_inputs": ["fasta"],
            "assets": {
                "kallisto_index": "kallisto_index"
                },
            "command_list": [
                "kallisto index -i {fasta} {asset_outfolder}/{genome}_kallisto_index.idx"
                ] 
        },
        "gtf_anno": {
            "required_inputs": ["gtf"],
            "assets": {
                "gtf_anno": "gtf_anno"
                },
            "command_list": [
                "cp {gtf} {asset_outfolder}/{genome}.gtf",
                ] 
        },
        "epilog_index": {
            "required_inputs": ["fasta", "context"],
            "assets": {
                "epilog_index": "epilog_index"
                },
            "command_list": [
                "epilog index -i {fasta} -o {asset_outfolder}/{genome}_{context}.tsv -s {context} -t"
                ] 
        }
    }


    def build_asset(genome, asset_key, asset_build_package, specific_args):
        _LOGGER.debug("Asset build package: " + str(asset_build_package))

        asset_outfolder = os.path.join(outfolder, asset_key)
        asset_vars = {"genome": genome,
                      "asset": asset_key,
                      "asset_outfolder": asset_outfolder}
        asset_vars.update(specific_args)


        print(str([x.format(**asset_vars) for x in asset_build_package["command_list"]]))

        tk.make_dir(asset_outfolder)
        target = os.path.join(asset_outfolder, "build_complete.flag")
        command_list_populated = [x.format(**asset_vars) for x in asset_build_package["command_list"]]

        touch_target = "touch {target}".format(target=target)
        command_list_populated.append(touch_target)
        
        _LOGGER.debug("Command list populated: " + str(command_list_populated))

        pm.run(command_list_populated, target, container=pm.container)
        # Add index information to rgc
        for asset_key, relative_path in asset_build_package["assets"].items():
            rgc.update_genomes(genome, asset_key, {"path": relative_path.format(**asset_vars)})

        # Write the updated refgenie genome configuration
        rgc.write()


    pm = pypiper.PipelineManager(name="refgenie", outfolder=outfolder, args=args)
    tk = pypiper.NGSTk(pm=pm)
    tools = pm.config.tools  # Convenience alias
    index = pm.config.index
    param = pm.config.param

    container = None
    if args.docker:
        # Set up some docker stuff
        pm.get_container("nsheff/refgenie", outfolder)

    for asset_key in args.asset:
        if asset_key in asset_build_packages.keys():
            asset_build_package = asset_build_packages[asset_key]
            _LOGGER.debug(specific_args)
            for required_input in asset_build_package["required_inputs"]:
                if not specific_args[required_input]:
                    raise ValueError("Argument '{}' is required to build asset '{}', but not provided".format(required_input, asset_key))
            build_asset(args.genome, asset_key, asset_build_package, specific_args)
        else:
            _LOGGER.warn("Recipe does not exist for asset '{}'".format(asset_key))


    if False:
        # pm.make_sure_path_exists(outfolder)
        conversions = {}
        conversions[".2bit"] = "twoBitToFa {INPUT} {OUTPUT}"
        conversions[".gz"] = tk.ziptool + " -cd {INPUT} > {OUTPUT}"

        # Copy fasta file to genome folder structure
        local_raw_fasta = genome + ".fa"
        raw_fasta = os.path.join(outfolder, local_raw_fasta)

        input_fasta, cmd = copy_or_download_file(args.fasta, outfolder)
        pm.run(cmd, input_fasta)

        cmd = convert_file(input_fasta, raw_fasta, conversions)
        if cmd:
            pm.run(cmd, raw_fasta, container=pm.container)


    # Copy annotation file (if any) to folder structure
    if args.gtf:
        annotation_file_unzipped = os.path.join(outfolder, genome + ".gtf")
        annotation_file, cmd = copy_or_download_file(args.gtf, outfolder)
        pm.run(cmd, annotation_file)

        cmd = convert_file(annotation_file, annotation_file_unzipped, conversions)
        pm.run(cmd, annotation_file_unzipped)

    #   cmd = "cp " + args.gtf + " " + annotation_file
    #   cmd2 = tk.ziptool + " -d " + annotation_file
    #   pm.run([cmd, cmd2], annotation_file_unzipped)

    else:
        _LOGGER.debug("* No GTF gene annotations provided. Skipping this step.")


    # # Bowtie indexes
    # if index.bowtie2:
    #     asset_key = "indexed_bowtie2"
    #     folder = os.path.join(outfolder, asset_key)
    #     tk.make_dir(folder)
    #     target = os.path.join(folder, "completed.flag")
    #     cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
    #     cmd2 = tools.bowtie2build + " " + raw_fasta + " " + os.path.join(folder, genome)
    #     cmd3 = "touch " + target
    #     pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
    #     # Add index information to rgc
    #     rgc.update_genomes(genome, asset_key, path_data(folder, rgc))

    #     # Write the updated refgenie genome configuration
    #     rgc.write()


    # # Bismark index - bowtie2
    # if index.bismark_bt2:
    #     asset_key = "indexed_bismark_bt2"
    #     folder = os.path.join(outfolder, asset_key)
    #     tk.make_dir(folder)
    #     target = os.path.join(folder, "completed.flag")
    #     cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
    #     cmd2 = tools.bismark_genome_preparation + " --bowtie2 " + folder
    #     cmd3 = "touch " + target
    #     pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
    #     rgc.update_genomes(genome, asset_key, path_data(folder, rgc))
    #     rgc.write()

    # # Bismark index - bowtie1
    # if index.bismark_bt1:
    #     asset_key = "indexed_bismark_bt1"
    #     folder = os.path.join(outfolder, asset_key)
    #     tk.make_dir(folder)
    #     target = os.path.join(folder, "completed.flag")
    #     cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
    #     cmd2 = tools.bismark_genome_preparation + " " + folder
    #     cmd3 = "touch " + target
    #     pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
    #     rgc.update_genomes(genome, asset_key, path_data(folder, rgc))
    #     rgc.write()

    # # Epilog meth calling
    # if index.epilog:
    #     asset_key = "indexed_epilog"
    #     folder = os.path.join(outfolder, asset_key)
    #     tk.make_dir(folder)
    #     target = os.path.join(folder, "completed.flag")
    #     cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
    #     cmd2 = tools.epilog_indexer + " -i " + raw_fasta
    #     cmd2 += " -o " + os.path.join(folder, genome + "_" + param.epilog.context + ".tsv")
    #     cmd2 += " -s " + param.epilog.context  # context
    #     cmd2 += " -t"
    #     cmd3 = "touch " + target
    #     pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
    #     rgc.update_genomes(genome, asset_key, path_data(folder, rgc))
    #     rgc.write()

    # if index.hisat2:
    #     asset_key = "indexed_hisat2"
    #     folder = os.path.join(outfolder, asset_key)
    #     tk.make_dir(folder)
    #     target = os.path.join(folder, "completed.flag")
    #     cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
    #     cmd2 = tools.hisat2build + " " + raw_fasta + " " + os.path.join(folder, genome)
    #     cmd3 = "touch " + target
    #     pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
    #     rgc.update_genomes(genome, asset_key, path_data(folder, rgc))
    #     rgc.write()

    # # Kallisto should index transcriptome
    # # So it doesn't make sense to run these at the same time as the others.
    # if index.kallisto:
    #     asset_key = "indexed_kallisto"
    #     folder = os.path.join(outfolder, asset_key)
    #     tk.make_dir(folder)
    #     target = os.path.join(folder, "completed.flag")
    #     cmd2 = tools.kallisto + " index -i " + os.path.join(folder, genome + "_kallisto_index.idx")
    #     cmd2 += " " + raw_fasta
    #     cmd3 = "touch " + target
    #     pm.run([cmd2, cmd3], target, container=pm.container)
    #     rgc.update_genomes(genome, asset_key, path_data(folder, rgc))
    #     rgc.write()

    pm.stop_pipeline()


def refgenie_init(genome_config_path, genome_server=DEFAULT_SERVER):
    """
    Initialize a genome config file.
    
    :param str genome_config_path: path to genome configuration file to 
        create/initialize
    :param st genome_server: URL for a server
    """

    # Set up default 
    rgc = RefGenConf(OrderedDict({
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


def _exec_list(rgc, remote):
    if remote:
        pfx = "Remote"
        assemblies, assets = rgc.list_remote()
    else:
        pfx = "Local"
        assemblies, assets = rgc.list_local()
    return pfx, assemblies, assets


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

    _LOGGER.debug("Args: {}".format(args))

    if not args.command:
        parser.print_help()
        _LOGGER.error("No command given")
        sys.exit(1)

    gencfg = yacman.select_config(
        args.genome_config, CFG_ENV_VARS,
        check_exist=not args.command == INIT_CMD,  on_missing=lambda fp: fp)
    if gencfg is None:
        raise MissingGenomeConfigError(args.genome_config)
    _LOGGER.debug("Determined genome config: {}".format(gencfg))

    if args.command == INIT_CMD:
        _LOGGER.info("Initializing refgenie genome configuration")
        _writeable(os.path.dirname(gencfg), strict_exists=True)
        refgenie_init(gencfg, args.genome_server)
        sys.exit(0)

    rgc = RefGenConf(gencfg)

    if args.command == BUILD_CMD:
        refgenie_build(rgc, args)

    elif args.command == GET_ASSET_CMD:
        _LOGGER.debug("getting asset: '{}/{}'".format(args.genome, args.asset))
        print(" ".join([rgc.get_asset(args.genome, asset) for asset in args.asset]))
        return

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
        rgc.pull_asset(args.genome, args.asset, gencfg, unpack=not args.no_untar)

    elif args.command in [LIST_LOCAL_CMD, LIST_REMOTE_CMD]:
        pfx, genomes, assets = _exec_list(rgc, args.command == LIST_REMOTE_CMD)
        _LOGGER.info("{} genomes: {}".format(pfx, genomes))
        _LOGGER.info("{} assets:\n{}".format(pfx, assets))


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
