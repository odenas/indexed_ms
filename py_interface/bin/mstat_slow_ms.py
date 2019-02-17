#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import sys
import argparse
from typing import Iterator

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def iter_prefixes(t: str, i: int) -> Iterator[str]:
    """
    Prefixes of t[i:], from longest to shortest
    """

    for j in reversed(range(i+1, len(t) + 1)):
        yield t[i:j]


def ms(t: str, s: str, i: int) -> int:
    """
    the longest prefix of t[i:] that occurs in s
    """

    for prefix in iter_prefixes(t, i):
        if prefix in s:
            return len(prefix)
    return 0


def main(opt):
    with open(opt.t_path) as fd:
        t = fd.read().rstrip()
    with open(opt.s_path) as fd:
        s = fd.read().rstrip()

    print(" ".join([str(ms(t, s, i)) for i in range(len(t))]),  "")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument('s_path', type=str,
                            help='Path of the S string i.e., index')
    arg_parser.add_argument('t_path', type=str,
                            help='Path of the T string i.e., query.')
    sys.exit(main(arg_parser.parse_args()))
