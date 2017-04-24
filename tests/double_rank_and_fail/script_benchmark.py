#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs while recording space/time usage
"""

import logging
import sys
import argparse
import os

from utils import MsInterface, verbose_args, get_output

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    base_dict = dict(lazy_wl=opt.lazy_wl,
                     load_cst=opt.load_cst,
                     space_usage=True, time_usage=True,
                     answer=False,
                     runs_progress=opt.runs_progress,
                     ms_progress=opt.ms_progress)
    for i, pref in enumerate(opt.prefixes):
        bpath = os.path.join(opt.base_dir, pref)
        command = MsInterface.ms_command_from_dict(dict(s_path=bpath + ".s",
                                                        t_path=bpath + ".t",
                                                        **base_dict))

        for j in range(opt.repeat):
            with open(opt.output, 'a') as fd:
                res = get_output(command)
                if i + j == 0:
                    fd.write(res[0] + ",label,b_path\n")
                for line in res[1:]:
                    fd.write(line.replace(" ", "") +
                             ("," + opt.label) +
                             ("," + os.path.basename(bpath)) +
                             "\n")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")

    for k in ('lazy_wl', 'rank_fail', 'load_cst'):
        args, kwargs = MsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)

    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    arg_parser.add_argument("--label", type=str, default='default',
                            help="ad a label col to output")
    arg_parser.add_argument("--repeat", type=int, default=1,
                            help="repeat each experiment")
    verbose_args(arg_parser)
    sys.exit(main(arg_parser.parse_args()))
