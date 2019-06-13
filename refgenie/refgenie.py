#!/usr/bin/env python

from argparse import ArgumentParser
from collections import OrderedDict
import os
import re
import sys

from ._version import __version__
from .exceptions import MissingGenomeConfigError, MissingFolderError

import logmuse
import pypiper
import refgenconf
from refgenconf import select_genome_config, RefGenConf
from refgenconf.const import *
from ubiquerg import is_url

_LOGGER = None

BUILD_CMD = "build"
INIT_CMD = "init"
PULL_CMD = "pull"
LIST_LOCAL_CMD = "list"
LIST_REMOTE_CMD = "listr"
GET_ASSET_CMD = "seek"


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
        sps[cmd].add_argument(
            '-c', '--genome-config', dest="genome_config",
            help="Path to local genome configuration file, to read from and/or "
                 "to create or update, depending on the operation")

    sps[INIT_CMD].add_argument('-s', '--genome-server', default=DEFAULT_SERVER,
                help="URL to use for the genome_server attribute in config file."
                " Defaults : {}".format(DEFAULT_SERVER))
    sps[BUILD_CMD] = pypiper.add_pypiper_args(
        sps[BUILD_CMD], groups=None, args=["recover", "config"])

    # Add any arguments specific to subcommands.

    sps[BUILD_CMD].add_argument(
        '-i', '--input', dest='input', required=True,
        help='Local path or URL to genome sequence file in .fa, .fa.gz, '
             'or .2bit format.')
    sps[BUILD_CMD].add_argument(
        '-n', '--name', dest='name', required=False,
        help='Name of the genome to build. If omitted, refgenie will use the '
             'basename of the file specified in --input.')
    sps[BUILD_CMD].add_argument(
        '-a', '--annotation', dest='annotation', required=False,
        help='Path to GTF gene annotation file.')
    sps[BUILD_CMD].add_argument(
        "-d", "--docker", action="store_true",
        help="Run all commands in the refgenie docker container.")
    sps[BUILD_CMD].add_argument(
        '-o', '--outfolder', dest='outfolder', required=False, default=None,
        help='Override the default path to genomes folder, which is to use the '
             'genome_folder attribute in the genome configuration file.')

    for cmd in [PULL_CMD, GET_ASSET_CMD]:
        sps[cmd].add_argument(
            "-g", "--genome", required=True,
            help="Reference assembly ID, e.g. mm10")
        sps[cmd].add_argument(
            "-a", "--asset", required=True, nargs='+',
            help="Name of asset, a key in a genome config file")

    sps[PULL_CMD].add_argument(
        "-u", "--no-untar", action="store_true",
        help="Do not extract tarballs.")

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


def convert_file(input_file, output_file, conversions):
    """
    Given an input file, output file, and a list of conversions, gives the appropriate output file.
    
    :param str output_file: Path to local output file you want to create
    :param dict conversions: A dictionary of shell commands to convert files of a given type.
    """
    form = {"INPUT": input_file, "OUTPUT": output_file}
    _, ext = os.path.splitext(input_file)
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

    if args.name:
        genome_name = args.name
    else:
        genome_name = os.path.basename(args.input)
        # eliminate extensions to get canonical genome name.
        for strike in [".fasta.gz$", ".fa.gz$", ".fasta$", ".fa$", ".gz$", ".2bit$"]:
            genome_name = re.sub(strike, "", genome_name)

    _LOGGER.info("Using genome name: {}".format(genome_name))

    if not hasattr(args, "outfolder") or not args.outfolder:
        # Default to genome_folder
        _LOGGER.debug("No outfolder provided, using genome config.")
        args.outfolder = rgc.genome_folder

    outfolder = os.path.abspath(os.path.join(args.outfolder, genome_name))
    if not _writeable(outfolder):
        _LOGGER.error("Insufficient permissions to write to output folder: {}".
                      format(outfolder))
        return

    _LOGGER.info("Output to: {} {} {}".format(genome_name, args.outfolder, outfolder))
    _LOGGER.debug("Default config file: {}".format(default_config_file()))

    if args.config_file and not os.path.isfile(args.config_file):
        _LOGGER.debug("Config file path isn't a file: {}".
                      format(args.config_file))
        args.config_file = default_config_file()

    pm = pypiper.PipelineManager(name="refgenie", outfolder=outfolder, args=args)
    tk = pypiper.NGSTk(pm=pm)
    tools = pm.config.tools  # Convenience alias
    index = pm.config.index
    param = pm.config.param

    # pm.make_sure_path_exists(outfolder)
    conversions = {}
    conversions[".2bit"] = "twoBitToFa {INPUT} {OUTPUT}"
    conversions[".gz"] = tk.ziptool + " -cd {INPUT} > {OUTPUT}"

    # Copy fasta file to genome folder structure
    local_raw_fasta = genome_name + ".fa"
    raw_fasta = os.path.join(outfolder, local_raw_fasta)

    input_file, cmd = copy_or_download_file(args.input, outfolder)
    pm.run(cmd, input_file)

    container = None
    if args.docker:
        # Set up some docker stuff
        pm.get_container("nsheff/refgenie", outfolder)
        # of = os.path.abspath(outfolder)
        # outfolder = of
        # cmd = "docker run -itd"
        # cmd += " -v " + of + ":" + of
        # cmd += " nsheff/refgenie"
        # container = pm.checkprint(cmd).rstrip()
        # _LOGGER.info("Using docker container: " + container)
        # pm.atexit_register(remove_container, container)

    cmd = convert_file(input_file, raw_fasta, conversions)
    if cmd:
        pm.run(cmd, raw_fasta, container=pm.container)

    cmd = tools.samtools + " faidx " + raw_fasta
    pm.run(cmd, raw_fasta + ".fai", container=pm.container)

    # Determine chromosome sizes
    fai_file = raw_fasta + ".fai"
    # symlinks should be relative so folders are portable.
    loc_chrom_sizes_file = genome_name + ".chrom.sizes"
    chrom_sizes_file = os.path.join(outfolder, loc_chrom_sizes_file)
    chrom_sizes_alias = os.path.join(outfolder, genome_name + ".chromSizes")
    cmd = "cut -f 1,2 " + fai_file + " > " + chrom_sizes_file
    cmd2 = "ln -s " + loc_chrom_sizes_file + " " + chrom_sizes_alias
    pm.run([cmd, cmd2], chrom_sizes_alias, container=pm.container)

    # Copy annotation file (if any) to folder structure
    if args.annotation:
        annotation_file_unzipped = os.path.join(outfolder, genome_name + ".gtf")
        annotation_file, cmd = copy_or_download_file(args.annotation, outfolder)
        pm.run(cmd, annotation_file)

        cmd = convert_file(annotation_file, annotation_file_unzipped, conversions)
        pm.run(cmd, annotation_file_unzipped)

    #   cmd = "cp " + args.annotation + " " + annotation_file
    #   cmd2 = tk.ziptool + " -d " + annotation_file
    #   pm.run([cmd, cmd2], annotation_file_unzipped)

    else:
        _LOGGER.debug("* No GTF gene annotations provided. Skipping this step.")

    def path_data(root, c):
        return {"path": os.path.relpath(root, c.genome_folder)}

    # Bowtie indexes
    if index.bowtie2:
        asset_key = "indexed_bowtie2"
        folder = os.path.join(outfolder, asset_key)
        tk.make_dir(folder)
        target = os.path.join(folder, "completed.flag")
        cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
        cmd2 = tools.bowtie2build + " " + raw_fasta + " " + os.path.join(folder, genome_name)
        cmd3 = "touch " + target
        pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
        # Add index information to rgc
        rgc.update_genomes(genome_name, asset_key, path_data(folder, rgc))

        # Write the updated refgenie genome configuration
        rgc.write()

    # Bismark index - bowtie2
    if index.bismark_bt2:
        asset_key = "indexed_bismark_bt2"
        folder = os.path.join(outfolder, asset_key)
        tk.make_dir(folder)
        target = os.path.join(folder, "completed.flag")
        cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
        cmd2 = tools.bismark_genome_preparation + " --bowtie2 " + folder
        cmd3 = "touch " + target
        pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
        rgc.update_genomes(genome_name, asset_key, path_data(folder, rgc))
        rgc.write()

    # Bismark index - bowtie1
    if index.bismark_bt1:
        asset_key = "indexed_bismark_bt1"
        folder = os.path.join(outfolder, asset_key)
        tk.make_dir(folder)
        target = os.path.join(folder, "completed.flag")
        cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
        cmd2 = tools.bismark_genome_preparation + " " + folder
        cmd3 = "touch " + target
        pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
        rgc.update_genomes(genome_name, asset_key, path_data(folder, rgc))
        rgc.write()

    # Epilog meth calling
    if index.epilog:
        asset_key = "indexed_epilog"
        folder = os.path.join(outfolder, asset_key)
        tk.make_dir(folder)
        target = os.path.join(folder, "completed.flag")
        cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
        cmd2 = tools.epilog_indexer + " -i " + raw_fasta
        cmd2 += " -o " + os.path.join(folder, genome_name + "_" + param.epilog.context + ".tsv")
        cmd2 += " -s " + param.epilog.context  # context
        cmd2 += " -t"
        cmd3 = "touch " + target
        pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
        rgc.update_genomes(genome_name, asset_key, path_data(folder, rgc))
        rgc.write()

    if index.hisat2:
        asset_key = "indexed_hisat2"
        folder = os.path.join(outfolder, asset_key)
        tk.make_dir(folder)
        target = os.path.join(folder, "completed.flag")
        cmd1 = "ln -sf ../" + local_raw_fasta + " " + folder
        cmd2 = tools.hisat2build + " " + raw_fasta + " " + os.path.join(folder, genome_name)
        cmd3 = "touch " + target
        pm.run([cmd1, cmd2, cmd3], target, container=pm.container)
        rgc.update_genomes(genome_name, asset_key, path_data(folder, rgc))
        rgc.write()

    # Kallisto should index transcriptome
    # So it doesn't make sense to run these at the same time as the others.
    if index.kallisto:
        asset_key = "indexed_kallisto"
        folder = os.path.join(outfolder, asset_key)
        tk.make_dir(folder)
        target = os.path.join(folder, "completed.flag")
        cmd2 = tools.kallisto + " index -i " + os.path.join(folder, genome_name + "_kallisto_index.idx")
        cmd2 += " " + raw_fasta
        cmd3 = "touch " + target
        pm.run([cmd2, cmd3], target, container=pm.container)
        rgc.update_genomes(genome_name, asset_key, path_data(folder, rgc))
        rgc.write()

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
        assemblies, assets = rgc.genomes_str(), rgc.assets_str()
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

    gencfg = select_genome_config(args.genome_config)
    if gencfg is None:
        raise MissingGenomeConfigError(args.genome_config)
    _LOGGER.info("Determined genome config: {}".format(gencfg))

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
        return " ".join([rgc.get_asset(args.genome, asset) for asset in args.asset])

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
