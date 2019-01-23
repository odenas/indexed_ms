#!/usr/bin/env python


"""
compute matching statistics
"""

import logging
import sys
import argparse

from mstat.text_stats import block_sample_iter


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def main(opt):
    with open(opt.input) as fd:
        txt = fd.read().strip()
        for i, (idx, sample) in enumerate(block_sample_iter(txt, opt.block_size, opt.n)):
            print(sample)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("input", type=str, help="input")
    arg_parser.add_argument("block_size", type=int, help="block size")
    arg_parser.add_argument("n", type=int, help="number of samples")
    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    sys.exit(main(arg_parser.parse_args()))

