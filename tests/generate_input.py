#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import sys
import argparse
import os

from utils import create_ms_input


LG = logging.getLogger(__name__)


def main(opt):
    create_ms_input(opt.input_type,
                 (opt.t_path, opt.t_len), (opt.s_path, opt.s_len),
                 opt.source, opt.mutation_period)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("input_type", choices=('random', 'file', 'mutation'),
                            help="how to generate.")
    arg_parser.add_argument("source", type=str,
                            help="alphabet for random/mutation, in_file else")

    arg_parser.add_argument('s_path', type=str, help="The S string")
    arg_parser.add_argument("s_len", type=int, help="length of t")
    arg_parser.add_argument('t_path', type=str, help="The T string i.e. the query")
    arg_parser.add_argument("t_len", type=int, help="length of s")

    arg_parser.add_argument("--mutation_period", type=int, default=5,
                            help="mutation period")
    sys.exit(main(arg_parser.parse_args()))
