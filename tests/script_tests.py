#!/usr/bin/env python
"""
Test whether programss produce the same output
"""

import subprocess
import logging
import sys
from difflib import ndiff
import argparse

from generate_strings import Input

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")


def main(opt):
    for i in range(opt.n):
        LG.info("testing on %d", i)
        I = Input._parse(str(i), opt.base_dir)
        res1 = get_output("python " + opt.p1 + " " + " ".join(I))
        res2 = get_output(opt.p2 + " " + " ".join(I))

        for i, (o1, o2) in enumerate(zip(res1, res2)):
            for line in ndiff(o1.split(), o2.split()):
                if line[:2] == '  ':
                    print "OK"
                else:
                    print line
            print


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("p1", type=str, help="program 1")
    arg_parser.add_argument("p2", type=str, help="program 2")
    arg_parser.add_argument("n", type=int, help="input prefixes 0..n")
    arg_parser.add_argument("--base_dir", type=str, help="base dir")
    sys.exit(main(arg_parser.parse_args()))
