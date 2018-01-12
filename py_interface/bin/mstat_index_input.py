#!/usr/bin/env python

"""
dump a string into a text file
"""


import sys
import logging
from collections import Counter
import random

import numpy as np

from mstat.dataset import rnd_textfile, rep_textfile
from mstat.interface import default_arg_parser, get_output, DumpCstInterface, DumpMaxrepInterface


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def main(opt):
    seed = (random.randint(0, 10000000) if opt.seed is None else opt.seed)
    LG.warning("SEED: %d", seed)
    random.seed(seed)
    np.random.seed(seed)

    if opt.itype == 'rnd':
        rnd_textfile(opt.output, opt.length, Counter(opt.alp))
    else:
        rep_textfile(opt.output,
                     opt.length, opt.num_blocks, opt.num_mut,
                     Counter(opt.alp))
    if not opt.skip_cst:
        get_output(DumpCstInterface.command_from_dict({'s_path': opt.output}))
    if not opt.skip_maxrep:
        get_output(DumpMaxrepInterface.command_from_dict({'s_path': opt.output}))


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)
    arg_parser.add_argument("itype", choices=('rnd', 'rep'),
                            help="Type of dataset")
    arg_parser.add_argument("length", type=int, help="string length.")
    arg_parser.add_argument('alp', type=str,
                            help="alphabet. repeat symbols to increase their frequency")
    arg_parser.add_argument('output', type=str, help='Outfle')
    arg_parser.add_argument('--skip_cst', action='store_true',
                            default=False, help='dont generate cst')
    arg_parser.add_argument('--skip_maxrep', action='store_true',
                            default=False, help='dont generate the maxrep vector')
    arg_parser.add_argument("--num_blocks", type=int,
                            help="nr. of blocks for repeat str.")
    arg_parser.add_argument("--num_mut", type=int,
                            help="nr. of mutations for repeat str.")
    sys.exit(main(arg_parser.parse_args()))
