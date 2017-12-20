#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import os
import sys
import argparse

sys.path.append("../../py_interface")
from bin_interfaces import default_arg_parser, CommandLineArguments, CommonArguments


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)

class IStats(CommandLineArguments):
    EXEC_PATH = "../../fast_ms/bin/input_stats.x"


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)
    LG.debug("running on S=%s, T=%s", opt.s_path, opt.t_path)
    command = InputStatsInterface.command_from_dict(vars(opt))
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

    for k in InputStatsInterface.params:
        args, kwargs = InputStatsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    arg_parser.add_argument("--output", type=str, default='/dev/stdout')
    verbose_args(arg_parser)

    sys.exit(main(arg_parser.parse_args()))
