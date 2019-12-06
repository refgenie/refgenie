#!/usr/bin/env python
"""
Each iGenome has the following nested directory structure:
    Species/
    Source/
    Build/
    Annotation/ Sequence/
"""
from .refgenie import refgenie_add
from .exceptions import MissingGenomeConfigError

from ubiquerg import untar, mkabs
from yacman import select_config

import refgenconf
from glob import glob

import os
import argparse
import sys
import tarfile
from shutil import move


def build_argparser():
    """
    Build a parser for this tool

    :return argparse.ArgumentParser: constructed parser
    """
    parser = argparse.ArgumentParser(description='Integrates every asset from the downloaded iGenomes'
                                                 ' tarball/directory with Refgenie asset management system')
    parser.add_argument('-p', '--path', dest="path", type=str,
                        help='path to the desired genome tarball or directory to integrate', required=True)
    parser.add_argument('-g', '--genome', dest="genome", type=str,  help='name to be assigned to the selected genome',
                        required=True)
    parser.add_argument('-c', '--config', dest="config", type=str,
                        help="path to local genome configuration file. Optional if '{}' environment variable is set.".
                        format(", ".join(refgenconf.CFG_ENV_VARS)), required=False)
    return parser


def untar_or_copy(p, dest):
    """
    Depending on a kind of the provided path, either copy or extract it to the destination directory

    :param str p: path to the directory to be copied or tarball to be extracted
    :param str dest: where to extract file or copy dir
    :return bool: whether the process was successful
    """
    if os.path.exists(p):
        if os.path.isdir(p):
            fun = move
            dest = os.path.join(dest, os.path.basename(p))
        elif tarfile.is_tarfile(p):
            fun = untar
            print("Extracting '{}'".format(p))
        else:
            raise ValueError("Provided path is neither a directory nor a tar archive.")
        fun(p, dest)
        print("Moved '{}' to '{}'".format(p, dest))
        return True
    return False


def main():
    """ main workflow """
    parser = build_argparser()
    args, remaining_args = parser.parse_known_args()
    cfg = refgenconf.select_genome_config(filename=args.config, check_exist=True, strict_env=True)
    if not cfg:
        raise MissingGenomeConfigError(args.config)
    rgc = refgenconf.RefGenConf(filepath=cfg, writable=False)
    pths = [args.path, mkabs(args.path, rgc.genome_folder)]
    if not untar_or_copy(pths[0], os.path.join(rgc.genome_folder, args.genome)) \
            and not untar_or_copy(pths[1], os.path.join(rgc.genome_folder, args.genome)):
        raise OSError("Path '{}' does not exist. Tried: {}".format(args.path, " and ".join(pths)))
    path_components = [rgc.genome_folder] + [args.genome] + ["*"] * 3 + ["Sequence"]
    assets_paths = glob(os.path.join(*path_components))
    assert len(assets_paths) > 0, OSError("Your iGenomes directory is corrupted, more than one directory matched by {}."
                                          "\nMatched dirs: {}".format(os.path.join(*path_components),
                                                                      ", ".join(assets_paths)))
    assets_path = assets_paths[0]
    asset_names = [d for d in os.listdir(assets_path) if os.path.isdir(assets_path)]
    processed = []
    for a in asset_names:
        asset_dict = {"genome": args.genome, "asset": a, "tag": None, "seek_key": None}
        asset_path = os.path.relpath(os.path.join(assets_path, a), rgc.genome_folder)
        if refgenie_add(rgc, asset_dict, asset_path):
            processed.append("{}/{}".format(asset_dict["genome"], asset_dict["asset"]))
    print("Added assets: \n- {}".format("\n- ".join(processed)))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
