#!/usr/bin/env python


"""
compute matching statistics
"""


import logging
import sys
import argparse

from utils import MsCommand, get_output, InputSpec

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    ispec = InputSpec(opt.base_dir, opt.prefix)
    LG.info("running on %s", ispec)
    res = get_output(MsCommand.fast(ispec,
                                    opt.lazy_wl,
                                    space_usage=True, time_usage=True,
                                    answer=True,
                                    verb=opt.vv,
                                    path_to_exec=opt.fast_prg))
    print res[-1]


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument('base_dir', type=str, help="directory of input")
    arg_parser.add_argument('prefix', type=str, help="prefix")
    arg_parser.add_argument("--fast_prg", type=str, default='./fd_ms',
                            help="c++ program")
    arg_parser.add_argument("--lazy_wl", action='store_true',
                            default=False, help="get lazy winer links")
    arg_parser.add_argument("--v", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("--vv", action='store_true',
                            default=False, help="very verbose")
    sys.exit(main(arg_parser.parse_args()))
