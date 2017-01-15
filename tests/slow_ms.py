#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import sys
import argparse

from utils import InputSpec

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def iter_prefixes(t, i):
    "from longest to shortest"

    for j in reversed(range(i+1, len(t) + 1)):
        yield t[i:j]


def ms(t, s, i):
    """
    the longest prefix of t[i:] that occurs in s
    """
    for prefix in iter_prefixes(t, i):
        if prefix in s:
            return len(prefix)
    return 0


def main(opt):
    inp = InputSpec(opt.base_dir, opt.prefix)

    with open(inp.t_path) as fd:
        t = fd.read().rstrip()
    with open(inp.s_path) as fd:
        s = fd.read().rstrip()

    print opt.prefix, " ".join([str(ms(t, s, i)) for i in range(len(t))]),  ""


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument('base_dir', type=str, help="directory of input")
    arg_parser.add_argument('prefix', type=str, help="prefix")
    sys.exit(main(arg_parser.parse_args()))
