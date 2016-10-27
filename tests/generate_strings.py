#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import sys
import random
import argparse


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def get_string(n):
    return "".join([random.choice(['a', 'b', 'c']) for i in range(n)])


def main(opt):
    for pair in range(opt.n):
        len_t = random.randint(15, 20)
        len_s = random.randint(15, 20)

        t, s = get_string(len_t), get_string(len_s)
        assert set(s) == set(t)
        print t, s


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("n", type=int, help="how many string pairs")
    sys.exit(main(arg_parser.parse_args()))

