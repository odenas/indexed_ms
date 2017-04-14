#!/usr/bin/env python

"""
Test whether programs produce the same output
"""

import subprocess
import logging
import sys
import argparse
import os
from difflib import ndiff

from utils import MsInput, verbose_args, MsInterface

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def slow(exec_path, s_path, t_path):
    return ("python {exec_path} {s_path} {t_path}".format(**locals()))

def fast(exec_path, s_path, t_path, lazy_wl, nthr):
    lazy_wl_flag = ("--lazy_wl" if lazy_wl else "")
    return ("python {exec_path} {s_path} {t_path} --answer {lazy_wl_flag} --nthreads {nthr}"
            .format(**locals()))

def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")


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
        bpath = os.path.join(opt.base_dir, pref)
        s_path=bpath + ".s"
        t_path=bpath + ".t"

        res1 = get_output(fast(opt.fast_prg, s_path, t_path, opt.lazy_wl, opt.nthreads))
        res2 = get_output(slow(opt.slow_prg, s_path, t_path))

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
    arg_parser.add_argument("--fast_prg", type=str, default='fast_ms.py',
                            help="fast python program")

    for k in ('lazy_wl', 'nthreads'):
        args, kwargs = MsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    verbose_args(arg_parser)
    sys.exit(main(arg_parser.parse_args()))
