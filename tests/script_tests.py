#!/usr/bin/env python

"""
Test whether programs produce the same output
"""

import logging
import sys
import argparse
from difflib import ndiff

from utils import verbose_args, MsInterface, MsInput, get_output


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def slow(exec_path, ms_input):
    s_path, t_path = ms_input
    return ("python {exec_path} {s_path} {t_path}".format(**locals()))


def fast(opt, ms_input):
    params = dict(lazy_wl=opt.lazy_wl,
                  rank_fail=opt.rank_fail,
                  use_maxrep=opt.use_maxrep,
                  answer=True,
                  s_path=ms_input.s_path, t_path=ms_input.t_path)
    return MsInterface.ms_command_from_dict(params)


def check_res(res1, res2):
    error_lst = []
    for i, (o1, o2) in enumerate(zip(res1, res2)):
        for line in ndiff(o1.split(), o2.split()):
            if line[:2] == '  ':
                error_lst.append(False)
            else:
                error_lst.append(line)
    return error_lst


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    for pref in opt.prefixes:
        ms_input = MsInput.basedir_form(opt.base_dir, pref)
        LG.info("testing on input %s", ms_input)

        res1 = get_output(fast(opt, ms_input))
        res2 = get_output(slow(opt.slow_prg, ms_input))

        err_lst = check_res(res1, res2)
        if any(err_lst):
            for err in err_lst:
                if err:
                    print "\t" + err
        print


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")

    arg_parser.add_argument("--slow_prg", type=str, default='slow_ms.py',
                            help="slow python program")

    for k in ('lazy_wl', 'rank_fail', 'use_maxrep', 'nthreads'):
        args, kwargs = MsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    verbose_args(arg_parser)
    sys.exit(main(arg_parser.parse_args()))
