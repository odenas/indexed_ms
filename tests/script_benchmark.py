#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs while recording space/time usage
"""

import logging
import sys
import argparse

from utils import MsCommand, get_output, InputSpec, FDMS_PATH

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    for i, pref in enumerate(opt.prefixes):
        ispec = InputSpec(opt.base_dir, pref)
        LG.info("running on %s", ispec)
        command_str = MsCommand.fast(ispec,
                                     opt.lazy_wl, opt.sada,
                                     space_usage=True, time_usage=True,
                                     answer=False,
                                     verb=opt.vv,
                                     opt.runs_progress, opt.ms_progress,
                                     path_to_exec=opt.fast_prg)
        for j in range(opt.repeat):
            with open(opt.output, 'a') as fd:
                res = get_output(command_str)
                if i + j == 0:
                    fd.write(res[0] + ",label\n")
                for line in res[1:]:
                    fd.write(line.replace(" ", "") + ("," + opt.label) + "\n")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("--fast_prg", type=str, default=FDMS_PATH,
                            help="c++ program")
    arg_parser.add_argument("--lazy_wl", action='store_true',
                            default=False, help="get lazy winer links")
    arg_parser.add_argument("--sada", action='store_true',
                            default=False, help="Sadakane's stree")
    arg_parser.add_argument("--v", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("--vv", action='store_true',
                            default=False, help="very verbose")
    arg_parser.add_argument("--runs_progress", type=int, default=0,
                            help="nr. of progress msgs for runs construction")
    arg_parser.add_argument("--ms_progress", type=int, default=0,
                            help="nr. of progress msgs for ms construction")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")
    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    arg_parser.add_argument("--label", type=str, default='default',
                            help="ad a label col to output")
    arg_parser.add_argument("--repeat", type=int, default=1,
                            help="repeat each experiment")
    sys.exit(main(arg_parser.parse_args()))
