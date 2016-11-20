#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs
"""

import subprocess
import logging
import sys
from difflib import ndiff
import argparse

from generate_input import Input

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    for i, pref in enumerate(opt.prefixes):
        LG.info("testing on %s", pref)
        input_args = opt.base_dir + " " + pref
        res = get_output(opt.fast_prg + " " + input_args)

        for line in res[0 if i == 0 else 1:]:
            print line


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("--fast_prg", type=str, default='./fd_ms',
                            help="c++ program")
    arg_parser.add_argument("--v", action='store_true', default=False, help="verbose")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")
    sys.exit(main(arg_parser.parse_args()))
