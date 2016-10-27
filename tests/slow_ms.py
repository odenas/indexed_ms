#!/usr/bin/env python


"""
compute matching statistics

/Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms
"""


import logging
import sys
import argparse


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)

def iter_prefixes(t, i):
    "from longest to shortest"

    for j in reversed(range(i+1, len(t) + 1)):
        yield t[i:j]


def ms(t, s, i):
    for prefix in iter_prefixes(t, i):
        if prefix in s:
            return len(prefix)
    return 0


def main(opt):
    with open(opt.infile) as fd:
        for line in fd:
            t, s = line.strip().split()
            print t, s, "".join([str(ms(t, s, i)) for i in range(len(t))])


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("infile", type=str, help="file of string pairs")
    sys.exit(main(arg_parser.parse_args()))
