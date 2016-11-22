#!/usr/bin/env python
"""
Test whether programss produce the same output
"""

import subprocess
import logging
import sys
from difflib import ndiff
import argparse

from utils import MsCommand

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")


def check_res(res1, res2):
    error_lst = []
    for i, (o1, o2) in enumerate(zip(res1, res2)):
        for line in ndiff(o1.split(), o2.split()):
            if line[:2] == '  ':
                error_lst.append(False)
            else:
                error_lst.append(line)
    return error_lst


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    for pref in opt.prefixes:
        LG.info("testing on %s", pref)
        res1 = get_output(MsCommand.fast(opt.base_dir, pref,
                                         False, True, False,
                                         opt.fast_prg))
        res2 = get_output(MsCommand.slow(opt.base_dir, pref, opt.slow_prg))

        err_lst = check_res(res1, res2)
        if any(err_lst):
            for err in err_lst:
                if err:
                    print "\t" + err
        print


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("--fast_prg", type=str, default='./fd_ms',
                            help="c++ program")
    arg_parser.add_argument("--slow_prg", type=str, default='slow_ms.py',
                            help="python program")
    arg_parser.add_argument("--v", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("--vv", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")
    sys.exit(main(arg_parser.parse_args()))
