#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import os
import sys
import argparse
import subprocess

from utils import MsInput, verbose_args, MsInterface
MsInterface.FDMS_PATH = ("/Users/denas/Library/Developer/Xcode/DerivedData/"
                         "fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/"
                         "Products/Debug/input_stats")


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def get_output(command, *args, **kwargs):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    return res.strip().split("\n")


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)
    LG.debug("running on S=%s, T=%s", opt.s_path, opt.t_path)
    command = MsInterface.command_from_dict(vars(opt))
    LG.debug("running: " + str(command))
    res = get_output(command)
    with open(opt.output, 'a') as fd:
        path = "b_path"
        for i, line in enumerate(res):
            if i:
                path = os.path.basename(opt.s_path)[:-2]
            fd.write(line + ("," + path) + "\n")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    for k in MsInterface.params:
        args, kwargs = MsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)

    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    arg_parser.add_argument("--label", type=str, default='default',
                            help="ad a label col to output")
    verbose_args(arg_parser)

    sys.exit(main(arg_parser.parse_args()))
