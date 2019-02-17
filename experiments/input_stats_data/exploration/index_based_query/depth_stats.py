#!/usr/bin/env python

"""
generate s and t strings
"""


import logging
import sys
import os

from mstat import get_output
from mstat.interface import default_arg_parser


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def main(opt):
    for i, path in enumerate(opt.input):
        if i == 0:
            print(",".join(('depth', 'node_depth', 'path', 'num_mut')))
        res = get_output(" ".join((opt.exe, path)))
        for line in res:
            b_path = os.path.basename(path)
            num_mut = int(b_path.split("_")[1][:-2])
            print(",".join((line, os.path.basename(path), str(num_mut))))


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)

    arg_parser.add_argument("input", type=str, nargs="+", help="Input files.")
    arg_parser.add_argument("--exe", type=str, default='./node_iterator.x',
                            help="Type of dataset.")
    sys.exit(main(arg_parser.parse_args()))
