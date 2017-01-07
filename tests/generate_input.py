#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import sys
import argparse
import os

from utils import InputSpec, create_input


LG = logging.getLogger(__name__)


def main(opt):
    prefix = (os.path.basename(opt.source) +
              ("%d_%d" % (opt.len_t, opt.len_t if opt.base_text == "mutation"
                                     else opt.len_s)))
    input_spec = InputSpec(opt.base_dir, prefix)
    create_input(input_spec, opt.len_t, opt.len_s,
                 opt.source, opt.base_text,
                 mutation_period=opt.mutation_period)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("base_text", choices=('random', 'file', 'mutation'),
                            help="a filename or 'random'")
    arg_parser.add_argument("source", type=str, help="alphabet for random, in_file for non-random")
    arg_parser.add_argument("--base_dir", type=str, default="./input_data",
                            help="base dir")
    arg_parser.add_argument("--len_t", type=int, default=100, help="length of t")
    arg_parser.add_argument("--len_s", type=int, default=5, help="length of s")
    arg_parser.add_argument("--mutation_period", type=int, default=5,
                            help="mutation period")
    sys.exit(main(arg_parser.parse_args()))
