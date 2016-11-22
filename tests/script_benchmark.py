#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs
"""

import subprocess
import logging
import sys
import argparse

from utils import MsCommand

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    with open(opt.output, 'w') as fd:
        for i, pref in enumerate(opt.prefixes):
            LG.info("running on %s", pref)
            res = get_output(MsCommand.fast(opt.base_dir,
                                            pref,
                                            True,
                                            False,
                                            opt.vv,
                                            opt.fast_prg))

            for line in res[(0 if i == 0 else 1):]:
                fd.write(line + "\n")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("--fast_prg", type=str, default='./fd_ms',
                            help="c++ program")
    arg_parser.add_argument("--v", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("--vv", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")
    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    sys.exit(main(arg_parser.parse_args()))
