"""
Test whether programss produce the same output
"""

import subprocess
import logging
import sys
from difflib import ndiff
import argparse

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def get_output(args):
    if args[0].endswith(".py"):
        command = 'python ' + " ".join(args)
    else:
        command = " ".join(args)
    LG.info("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.info("got: " + res)
    return res.strip().split("\n")


def main(opt):
    res1 = get_output([opt.p1, opt.in_file])
    res2 = get_output([opt.p2, opt.in_file])

    # res1 = ["baabcaabba cacabbcacc 1123213211", "cbaaacba bcbbcba 32111321", "bbabbc acabbbca 213321"]
    # res2 = ["baabcaabba cacabbcacc 1123543221", "cbaaacba bcbbcba 32111321", "bbabbc acabbbca 243321"]

    for o1, o2 in zip(res1, res2):
        for line in ndiff(o1.split(), o2.split()):
            print line.strip()
        print


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("p1", type=str, help="program 1")
    arg_parser.add_argument("p2", type=str, help="program 2")
    arg_parser.add_argument("in_file", type=str, help="input")
    sys.exit(main(arg_parser.parse_args()))
