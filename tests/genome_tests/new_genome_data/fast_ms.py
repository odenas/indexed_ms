#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import sys
import argparse
import subprocess

from utils import MsInput, verbose_args, MsInterface


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def get_output(command, *args, **kwargs):
    LG.debug("running: " + str(command))
    res = subprocess.check_call(command, *args, shell=True, **kwargs)
    return res


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)
    LG.debug("running on S=%s, T=%s", opt.s_path, opt.t_path)
    command = MsInterface.ms_command_from_dict(vars(opt))
    LG.debug("running: " + str(command))
    subprocess.check_call(command, shell=True)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    for k in MsInterface.params:
        args, kwargs = MsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)

    verbose_args(arg_parser)

    sys.exit(main(arg_parser.parse_args()))
