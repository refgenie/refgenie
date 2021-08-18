from argparse import HelpFormatter

import pypiper
from refgenconf import __version__ as rgc_version
from ubiquerg import VersionInHelpParser

from ._version import __version__
from .const import *


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
        epilog=additional_description,
    )

    subparsers = parser.add_subparsers(dest="command")

    def add_subparser(cmd, msg, subparsers):
        return subparsers.add_parser(
            cmd,
            description=msg,
            help=msg,
            formatter_class=lambda prog: HelpFormatter(
                prog, max_help_position=40, width=90
            ),
        )

    sps = {}
    for cmd, desc in SUBPARSER_MESSAGES.items():
        sps[cmd] = add_subparser(cmd, desc, subparsers)
        # alias is nested and alias subcommands require config path
        if cmd == ALIAS_CMD:
            continue
        # It's required for init
        sps[cmd].add_argument(
            "-c",
            "--genome-config",
            required=(cmd == INIT_CMD),
            dest="genome_config",
            metavar="C",
            help="Path to local genome configuration file. Optional if {} environment variable is set.".format(
                ", ".join(CFG_ENV_VARS)
            ),
        )
        sps[cmd].add_argument(
            "--skip-read-lock",
            required=False,
            action="store_true",
            help="Whether the config file should not be locked for reading",
        )

    # upgrade: upgrade config and alter file structure to the target version
    sps[UPGRADE_CMD].add_argument(
        "-v",
        "--target-version",
        required=True,
        metavar="V",
        help="Target config version for the upgrade.",
    )
    sps[UPGRADE_CMD].add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Do not prompt before action, approve upfront.",
    )

    sps[INIT_CMD].add_argument(
        "-s",
        "--genome-server",
        nargs="+",
        metavar="S",
        default=[DEFAULT_SERVER],
        help=f"URL(s) to use for the {CFG_SERVERS_KEY} attribute in config file. Default: {DEFAULT_SERVER}.",
    )
    sps[INIT_CMD].add_argument(
        "-f",
        "--genome-folder",
        metavar="F",
        help="Absolute path to parent folder refgenie-managed assets.",
    )
    sps[INIT_CMD].add_argument(
        "-r",
        "--recipe-folder",
        metavar="F",
        help="Absolute path to parent folder refgenie-managed recipes.",
    )
    sps[INIT_CMD].add_argument(
        "-l",
        "--asset-class-folder",
        metavar="F",
        help="Absolute path to parent folder refgenie-managed asset classes.",
    )
    sps[INIT_CMD].add_argument(
        "-a",
        "--genome-archive-folder",
        metavar="F",
        help="Absolute path to parent archive folder refgenie-managed assets; used by refgenieserver.",
    )
    sps[INIT_CMD].add_argument(
        "-b",
        "--genome-archive-config",
        metavar="C",
        help="Absolute path to desired archive config file; used by refgenieserver.",
    )
    sps[INIT_CMD].add_argument(
        "-u",
        "--remote-url-base",
        metavar="U",
        help="URL to use as an alternative, remote archive location; used by refgenieserver.",
    )
    sps[INIT_CMD].add_argument(
        "-j",
        "--settings-json",
        metavar="J",
        help="Absolute path to a JSON file with the key "
        "value pairs to inialize the configuration "
        "file with. Overwritten by itemized specifications.",
    )
    sps[BUILD_CMD] = pypiper.add_pypiper_args(
        sps[BUILD_CMD], groups=None, args=["recover", "config", "new-start"]
    )

    # Add any arguments specific to subcommands.

    sps[BUILD_CMD].add_argument(
        "--tag-description",
        required=False,
        default=None,
        type=str,
        help="Add tag level description (e.g. built with version 0.3.2).",
    )

    sps[BUILD_CMD].add_argument(
        "--genome-description",
        required=False,
        default=None,
        type=str,
        help="Add genome level description (e.g. The mouse mitochondrial genome, released in Dec 2013).",
    )

    sps[BUILD_CMD].add_argument(
        "-d",
        "--docker",
        action="store_true",
        help="Run all commands in the refgenie docker container.",
    )

    sps[BUILD_CMD].add_argument(
        "--pipeline-kwargs",
        nargs="+",
        action="append",
        required=False,
        default=None,
        help="Extra arguments to pass to the build pipeline. Format: arg_name=arg_val arg_name1=arg_val1",
    )

    sps[BUILD_CMD].add_argument(
        "--map",
        action="store_true",
        help="Run the map procedure: build assets and store the metadata in separate configs.",
    )

    sps[BUILD_CMD].add_argument(
        "--pull-parents",
        action="store_true",
        help="Automatically pull the default parent asset if required but not provided",
    )

    sps[BUILD_CMD].add_argument(
        "--preserve-map-configs",
        action="store_true",
        help="Do not remove the genome configuration files produced in the potential map step of building",
    )

    sps[BUILD_CMD].add_argument(
        "--reduce",
        action="store_true",
        help="Run the reduce procedure: gather the metadata produced with `refgenie build --map`.",
    )

    sps[BUILD_CMD].add_argument(
        "--assets",
        nargs="+",
        action="append",
        required=False,
        default=None,
        help="Override the default genome, asset and tag of the parents"
        " (e.g. fasta=hg38/fasta:default gtf=mm10/gencode_gtf:default).",
    )

    sps[BUILD_CMD].add_argument(
        "--files",
        nargs="+",
        action="append",
        required=False,
        default=None,
        help="Provide paths to the required files (e.g. fasta=/path/to/file.fa.gz).",
    )

    sps[BUILD_CMD].add_argument(
        "--params",
        nargs="+",
        action="append",
        required=False,
        default=None,
        help="Provide required parameter values (e.g. param1=value1).",
    )

    sps[BUILD_CMD].add_argument(
        "-v",
        "--volumes",
        nargs="+",
        required=False,
        default=None,
        help="If using docker, also mount these folders as volumes.",
    )

    sps[BUILD_CMD].add_argument(
        "-q",
        "--requirements",
        action="store_true",
        help="Show the build requirements for the specified asset and exit.",
    )

    sps[BUILD_CMD].add_argument(
        "--text",
        action="store_true",
        help="Can be used with '--requirements' to display the recipe requirements in text form.",
    )

    sps[BUILD_CMD].add_argument(
        "-r",
        "--recipe",
        required=False,
        default=None,
        type=str,
        help="Provide a recipe to use.",
    )

    recipe_subparser = sps[RECIPE_CMD]
    recipe_subsubparsers = recipe_subparser.add_subparsers(dest="subcommand")

    recipe_sps = {}
    for cmd, desc in RECIPE_SUBPARSER_MESSAGES.items():
        recipe_sps[cmd] = add_subparser(cmd, desc, recipe_subsubparsers)
        recipe_sps[cmd].add_argument(
            "-c",
            "--genome-config",
            required=False,
            dest="genome_config",
            metavar="C",
            help="Path to local genome configuration file. Optional if {} environment variable is set.".format(
                ", ".join(CFG_ENV_VARS)
            ),
        )
        recipe_sps[cmd].add_argument(
            "--skip-read-lock",
            required=False,
            action="store_true",
            help="Whether the config file should not be locked for reading",
        )

        if cmd not in [RECIPE_LIST_CMD, RECIPE_ADD_CMD, RECIPE_LIST_REMOTE_CMD]:
            recipe_sps[cmd].add_argument(
                "recipe_name",
                metavar="RECIPE_NAME",
                type=str,
                nargs=1,
                help="Recipe name to perform the action on.",
            )
        if cmd == RECIPE_ADD_CMD:
            recipe_sps[cmd].add_argument(
                "--source",
                required=True,
                type=str,
                help="Path/URL to the recipe to add.",
            )

        if cmd == RECIPE_PULL_CMD:
            recipe_sps[cmd].add_argument(
                "-f",
                "--force",
                required=False,
                action="store_true",
                help="Whether to force the pull even if the recipe exists.",
            )

            recipe_sps[cmd].add_argument(
                "--skip-asset-class",
                required=False,
                action="store_true",
                help="Whether to skip asset class pull.",
            )

        if cmd in [RECIPE_ADD_CMD, RECIPE_REMOVE_CMD]:
            recipe_sps[cmd].add_argument(
                "-f",
                "--force",
                required=False,
                action="store_true",
                help="Whether to force the action on the recipe.",
            )

        if cmd == RECIPE_SHOW_CMD:
            recipe_sps[cmd].add_argument(
                "-r",
                "--requirements",
                required=False,
                action="store_true",
                help="Whether to show the recipe requirements.",
            )
            recipe_sps[cmd].add_argument(
                "--remote",
                required=False,
                action="store_true",
                help="Whether to show the remote recipe.",
            )
            recipe_sps[cmd].add_argument(
                "-s",
                "--genome-server",
                nargs="+",
                required=False,
                metavar="S",
                help="One or more URLs to use. This information will not persist in the genome config file.",
            )
            recipe_sps[cmd].add_argument(
                "-p",
                "--append-server",
                action="store_true",
                help="Whether the provided servers should be appended to the list.",
            )

        if cmd == RECIPE_LIST_REMOTE_CMD:
            recipe_sps[cmd].add_argument(
                "-s",
                "--genome-server",
                nargs="+",
                required=False,
                metavar="S",
                help="One or more URLs to use. This information will not persist in the genome config file.",
            )
            recipe_sps[cmd].add_argument(
                "-p",
                "--append-server",
                action="store_true",
                help="Whether the provided servers should be appended to the list.",
            )

    asset_class_subparser = sps[ASSET_CLASS_CMD]
    asset_class_subsubparsers = asset_class_subparser.add_subparsers(dest="subcommand")

    asset_class_sps = {}
    for cmd, desc in ASSET_CLASS_SUBPARSER_MESSAGES.items():
        asset_class_sps[cmd] = add_subparser(cmd, desc, asset_class_subsubparsers)
        asset_class_sps[cmd].add_argument(
            "-c",
            "--genome-config",
            required=False,
            dest="genome_config",
            metavar="C",
            help="Path to local genome configuration file. Optional if {} environment variable is set.".format(
                ", ".join(CFG_ENV_VARS)
            ),
        )
        asset_class_sps[cmd].add_argument(
            "--skip-read-lock",
            required=False,
            action="store_true",
            help="Whether the config file should not be locked for reading",
        )
        if cmd not in [
            ASSET_CLASS_LIST_CMD,
            ASSET_CLASS_ADD_CMD,
            ASSET_CLASS_LIST_REMOTE_CMD,
        ]:
            asset_class_sps[cmd].add_argument(
                "asset_class_name",
                metavar="ASSET_CLASS_NAME",
                type=str,
                nargs=1,
                help="Asset class name to perform the action on.",
            )
        if cmd == ASSET_CLASS_ADD_CMD:
            asset_class_sps[cmd].add_argument(
                "--source",
                required=True,
                type=str,
                help="Path/URL to the asset class to add.",
            )

        if cmd in [ASSET_CLASS_ADD_CMD, ASSET_CLASS_REMOVE_CMD]:
            asset_class_sps[cmd].add_argument(
                "-f",
                "--force",
                required=False,
                action="store_true",
                help="Whether to force the action on the asset class.",
            )

        if cmd == ASSET_CLASS_PULL_CMD:
            asset_class_sps[cmd].add_argument(
                "-f",
                "--force",
                required=False,
                action="store_true",
                help="Whether to force the pull even if the asset class exists.",
            )

        if cmd == ASSET_CLASS_LIST_REMOTE_CMD:
            asset_class_sps[cmd].add_argument(
                "-s",
                "--genome-server",
                nargs="+",
                required=False,
                metavar="S",
                help="One or more URLs to use. This information will not persist in the genome config file.",
            )
            asset_class_sps[cmd].add_argument(
                "-p",
                "--append-server",
                action="store_true",
                help="Whether the provided servers should be appended to the list.",
            )
        if cmd == ASSET_CLASS_SHOW_CMD:
            asset_class_sps[cmd].add_argument(
                "-r",
                "--remote",
                required=False,
                action="store_true",
                help="Whether to show the remote asset class.",
            )
            asset_class_sps[cmd].add_argument(
                "-s",
                "--genome-server",
                nargs="+",
                required=False,
                metavar="S",
                help="One or more URLs to use. This information will not persist in the genome config file.",
            )
            asset_class_sps[cmd].add_argument(
                "-p",
                "--append-server",
                action="store_true",
                help="Whether the provided servers should be appended to the list.",
            )

    alias_subparser = sps[ALIAS_CMD]
    alias_subsubparsers = alias_subparser.add_subparsers(dest="subcommand")

    alias_sps = {}
    for cmd, desc in ALIAS_SUBPARSER_MESSAGES.items():
        alias_sps[cmd] = add_subparser(cmd, desc, alias_subsubparsers)
        alias_sps[cmd].add_argument(
            "-c",
            "--genome-config",
            required=False,
            dest="genome_config",
            metavar="C",
            help="Path to local genome configuration file. Optional if {} environment variable is set.".format(
                ", ".join(CFG_ENV_VARS)
            ),
        )
        alias_sps[cmd].add_argument(
            "--skip-read-lock",
            required=False,
            action="store_true",
            help="Whether the config file should not be locked for reading",
        )

    alias_sps[ALIAS_SET_CMD].add_argument(
        "-a",
        "--aliases",
        metavar="A",
        required=False,
        default=None,
        type=str,
        nargs="+",
        help="Aliases to set; single if the digest is to be retrieved from the server.",
    )
    alias_sps[ALIAS_SET_CMD].add_argument(
        "-d",
        "--digest",
        metavar="D",
        required=False,
        type=str,
        help="Digest to set; leave out if the digest is to be retrieved from the server.",
    )
    alias_sps[ALIAS_SET_CMD].add_argument(
        "-r",
        "--reset",
        action="store_true",
        help="Whether all the aliases should be removed prior to setting new ones.",
    )
    alias_sps[ALIAS_SET_CMD].add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Whether the action should be forced, if genome does not exist.",
    )

    alias_sps[ALIAS_REMOVE_CMD].add_argument(
        "-a",
        "--aliases",
        metavar="A",
        required=False,
        default=None,
        type=str,
        nargs="+",
        help="Aliases to remove.",
    )
    alias_sps[ALIAS_REMOVE_CMD].add_argument(
        "-d", "--digest", metavar="D", required=True, type=str, help="Digest to remove."
    )

    alias_sps[ALIAS_GET_CMD].add_argument(
        "-a",
        "--aliases",
        metavar="A",
        required=False,
        type=str,
        nargs="+",
        help="Aliases to get the digests for.",
    )

    sps[COMPARE_CMD].add_argument(
        "genome1",
        metavar="GENOME1",
        type=str,
        nargs=1,
        help="First genome for compatibility check.",
    )
    sps[COMPARE_CMD].add_argument(
        "genome2",
        metavar="GENOME2",
        type=str,
        nargs=1,
        help="Second genome for compatibility check.",
    )
    sps[COMPARE_CMD].add_argument(
        "-e",
        "--no-explanation",
        action="store_true",
        help="Do not print compatibility code explanation.",
    )
    sps[COMPARE_CMD].add_argument(
        "-f",
        "--flag-meanings",
        action="store_true",
        help="Display compatibility flag meanings.",
    )

    # add 'genome' argument to many commands
    for cmd in [
        PULL_CMD,
        GET_ASSET_CMD,
        GET_REMOTE_ASSET_CMD,
        BUILD_CMD,
        INSERT_CMD,
        REMOVE_CMD,
        GETSEQ_CMD,
        TAG_CMD,
        ID_CMD,
    ]:
        # genome is not required for listing actions
        sps[cmd].add_argument(
            "-g",
            "--genome",
            required=cmd in GETSEQ_CMD,
            metavar="G",
            help="Reference assembly ID, e.g. mm10.",
        )

    for cmd in LIST_REMOTE_CMD, LIST_LOCAL_CMD:
        sps[cmd].add_argument(
            "-g",
            "--genome",
            required=False,
            type=str,
            metavar="G",
            nargs="*",
            help="Reference assembly ID, e.g. mm10.",
        )

    for cmd in [
        PULL_CMD,
        GET_ASSET_CMD,
        GET_REMOTE_ASSET_CMD,
        BUILD_CMD,
        INSERT_CMD,
        REMOVE_CMD,
        TAG_CMD,
        ID_CMD,
    ]:
        build_arg_kwargs = dict(
            metavar="asset-registry-paths",
            type=str,
            nargs="+",
            help="One or more registry path strings that identify assets  (e.g. hg38/fasta or hg38/fasta:tag"
            + (
                " or hg38/fasta.fai:tag)."
                if cmd in [GET_ASSET_CMD, GET_REMOTE_ASSET_CMD]
                else ")."
            ),
        )
        # make asset-registry-path argument optional for build command
        # and require it manually in CLI when running a non-reduce build
        if cmd == BUILD_CMD:
            build_arg_kwargs.update({"nargs": "*", "default": None})
        sps[cmd].add_argument("asset_registry_paths", **build_arg_kwargs)

    sps[LIST_LOCAL_CMD].add_argument(
        "-r", "--recipes", action="store_true", help="List available recipes."
    )

    for cmd in [REMOVE_CMD, INSERT_CMD]:
        sps[cmd].add_argument(
            "-f",
            "--force",
            action="store_true",
            help="Do not prompt before action, approve upfront.",
        )

    sps[REMOVE_CMD].add_argument(
        "-a",
        "--aliases",
        action="store_true",
        help="Remove the genome alias if last asset for that genome is removed.",
    )
    force_group = sps[PULL_CMD].add_argument_group(
        title="Prompt handling",
        description="These flags configure the pull prompt responses.",
    )

    overwrite_group = force_group.add_mutually_exclusive_group()

    overwrite_group.add_argument(
        "--no-overwrite", action="store_true", help="Do not overwrite if asset exists."
    )

    overwrite_group.add_argument(
        "--force-overwrite", action="store_true", help="Overwrite if asset exists."
    )

    large_group = force_group.add_mutually_exclusive_group()

    large_group.add_argument(
        "--no-large", action="store_true", help="Do not pull archives over 5GB."
    )

    large_group.add_argument(
        "--pull-large",
        action="store_true",
        help="Pull any archive, regardless of its size.",
    )

    force_group.add_argument(
        "--size-cutoff",
        type=float,
        default=10,
        metavar="S",
        help="Maximum archive file size to download with no confirmation required (in GB, default: 10)",
    )

    force_group.add_argument(
        "-b",
        "--batch",
        action="store_true",
        help="Use batch mode: pull large archives, do no overwrite",
    )

    sps[INSERT_CMD].add_argument(
        "-p", "--path", required=True, metavar="P", help="Relative local path to asset."
    )

    sps[INSERT_CMD].add_argument(
        "-s",
        "--seek-keys",
        required=False,
        type=str,
        metavar="S",
        help="""
        String representation of a JSON object with seek_keys,
        e.g. '{"seek_key1": "file.txt"}'
        """,
    )

    sps[GETSEQ_CMD].add_argument(
        "-l",
        "--locus",
        required=True,
        help="Coordinates of desired sequence; e.g. 'chr1:50000-50200'.",
    )

    sps[GET_ASSET_CMD].add_argument(
        "-e",
        "--check-exists",
        required=False,
        action="store_true",
        help="Whether the returned asset path should be checked for existence on disk.",
    )

    sps[TAG_CMD].add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Do not prompt before action, approve upfront.",
    )

    group = sps[TAG_CMD].add_mutually_exclusive_group(required=True)

    group.add_argument("-t", "--tag", type=str, help="Tag to assign to an asset.")

    group.add_argument(
        "-d",
        "--default",
        action="store_true",
        help="Set the selected asset tag as the default one.",
    )

    sps[SUBSCRIBE_CMD].add_argument(
        "-r",
        "--reset",
        action="store_true",
        help="Overwrite the current list of server URLs.",
    )

    for cmd in [SUBSCRIBE_CMD, UNSUBSCRIBE_CMD]:
        sps[cmd].add_argument(
            "-s",
            "--genome-server",
            nargs="+",
            required=True,
            metavar="S",
            help="One or more URLs to {action} the {key} attribute in config file.".format(
                action="add to" if cmd == SUBSCRIBE_CMD else "remove from",
                key=CFG_SERVERS_KEY,
            ),
        )

    for cmd in [LIST_REMOTE_CMD, GET_REMOTE_ASSET_CMD, POPULATE_REMOTE_CMD]:
        sps[cmd].add_argument(
            "-s",
            "--genome-server",
            nargs="+",
            required=False,
            metavar="S",
            help="One or more URLs to use. "
            "This information will not persist in the genome config file.",
        )
        sps[cmd].add_argument(
            "-p",
            "--append-server",
            action="store_true",
            help="Whether the provided servers should be appended to the list.",
        )

    for cmd in [POPULATE_REMOTE_CMD, GET_REMOTE_ASSET_CMD]:
        sps[cmd].add_argument(
            "--remote-class",
            metavar="RC",
            type=str,
            default="http",
            help="Remote data provider class, e.g. 'http' or 's3'",
        )

    for cmd in [POPULATE_REMOTE_CMD, POPULATE_CMD]:
        sps[cmd].add_argument(
            "-f", "--file", metavar="F", help="File with registry paths to populate"
        )

    return parser
