#!/usr/bin/env python
"""
Each iGenomes has the following nested directory structure:
    Species/
    Source/
    Build/
    Annotation/ Sequence/
"""
from .refgenie import refgenie_add

from ubiquerg import untar, mkabs

import refgenconf
from glob import glob

import os
import argparse
import sys


def build_argparser():
    parser = argparse.ArgumentParser(description='Integrates every asset from the downloaded iGenomes tarball/directory '
                                                 'with Refgenie asset management system')
    parser.add_argument('-p', '--path', dest="path", type=str,  help='path to the desired genome tarball to integrate',
                        required=True)
    parser.add_argument('-g', '--genome', dest="genome", type=str,  help='name to be assigned to the selected genome',
                        required=True)
    parser.add_argument('-c', '--config', dest="config", type=str,  help='genome config', required=True)
    return parser


def main():
    parser = build_argparser()
    args, remaining_args = parser.parse_known_args()
    pths = [args.path, mkabs(args.path)]
    for p in pths:
        try:
            print("Extracting '{}' to '{}'".format(p, args.genome))
            untar(p, args.genome)
        except FileNotFoundError:
            continue
        else:
            break
        raise OSError("Provided path does not exist, tried: {}".format(" and".join(pths)))
    rgc = refgenconf.RefGenConf(args.config)
    assets_paths = glob(os.path.join(*([args.genome] + ["*"] * 3 + ["Sequence"])))
    assert len(assets_paths) == 1, OSError("your iGenomes directory is corrupted, more that one directory matched: {}".
                                           format(os.path.join(*([args.genome] + ["*"] * 3 + ["Sequence"]))))
    assets_path = assets_paths[0]
    asset_names = [dI for dI in os.listdir(assets_path) if os.path.isdir(assets_path)]
    for a in asset_names:
        asset_dict = {"genome": args.genome, "asset": a, "tag": None, "seek_key": None}
        asset_path = os.path.join(assets_path, a)
        print(asset_path)
        refgenie_add(rgc, asset_dict, asset_path)
    print("Asset added: {}/{}".format(asset_dict["genome"], asset_dict["asset"]))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
