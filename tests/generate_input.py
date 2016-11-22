#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import sys
import argparse
import os

from utils import Input


LG = logging.getLogger(__name__)


def main(opt):
    if opt.base_text == 'random':
        Input.random_build(opt)
    else:
        prefix = (os.path.basename(opt.in_file) +
                  ("%d_%d" % (opt.len_t, opt.len_s)))
        Input.file_build(opt.in_file, opt.len_t, opt.len_s,
                         prefix, opt.base_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("base_text", choices=('random', 'non_random'),
                            help="a filename or 'random'")
    arg_parser.add_argument("--base_dir", type=str, default="./input_data",
                            help="base dir")

    nrgrp = arg_parser.add_argument_group('options non-random data')
    nrgrp.add_argument("--in_file", type=str, help="base text")
    nrgrp.add_argument("--len_t", type=int, default=10,
                       help="t is the prefix of the input")
    nrgrp.add_argument("--len_s", type=int, default=5,
                       help="s follows t in the input")

    rgrp = arg_parser.add_argument_group('options random data')
    rgrp.add_argument("--n", type=int, help="how many string pairs")
    rgrp.add_argument("--l", type=int, default=7, help="length lbound")
    rgrp.add_argument("--L", type=int, default=10, help="length hbound")
    rgrp.add_argument("--alphabet", type=str, default='ab', help="alphabet")
    sys.exit(main(arg_parser.parse_args()))
