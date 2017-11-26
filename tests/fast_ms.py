#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import os
import sys
import argparse
import subprocess

#from utils import verbose_args, MsInterface, get_output
from bin_interfaces import default_arg_parser, FdMsInterface, get_output


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.verbose else logging.INFO)
    if opt.output != '/dev/stdout' and os.path.exists(opt.output):
        LG.error("output file (%s) exsts. Exiting ...", opt.output)
        return 1
    command = FdMsInterface.command_from_dict(vars(opt))

    pref = os.path.basename(opt.s_path).replace(".s", "")
    header_suff = ",".join(['label', opt.repeat_colname, 'b_path'])
    for j in range(opt.repeat):
        with open(opt.output, 'a') as fd:
            res = get_output(command)
            if j == 0:
                fd.write(res[0] + "," + header_suff + "\n")
            for line in res[1:]:
                fd.write(line.replace(" ", "") +
                         ("," + opt.label) +
                         ("," + str(j + 1)) +
                         ("," + pref) + "\n")


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)

    for k in FdMsInterface.params:
        args, kwargs = FdMsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    arg_parser.add_argument("--repeat", type=int, default=1,
                            help="repeat experiment")
    arg_parser.add_argument("--label", type=str, default='default',
                            help="ad a label column to output")
    arg_parser.add_argument("--repeat_colname", type=str, default='ntrial',
                            help="col label corresponding to the repeat nr.")
    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")

    sys.exit(main(arg_parser.parse_args()))
