#!/usr/bin/env python

"""
Test whether programs produce the same output
"""

import logging
import sys
import os

from mstat.interface import default_arg_parser, get_output, FdMsInterface
from mstat.dataset import InputPair


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def slow(load_dir, ms_input):
    s_path, t_path = ms_input
    if load_dir is None:
        return ("mstat_slow_ms.py {s_path} {t_path}".format(**locals()))
    load_path = os.path.join(load_dir, os.path.basename(s_path).replace(".s", ".mstat"))
    return ("cat {load_path}".format(**locals()))


def fast(opt, ms_input):
    params = dict(lazy_wl=opt.lazy_wl,
                  rank_fail=opt.rank_fail,
                  double_rank=opt.double_rank,
                  lca_parents=opt.lca_parents,
                  use_maxrep_vanilla=opt.use_maxrep_vanilla,
                  use_maxrep_rc=opt.use_maxrep_rc,
                  answer=opt.answer,
                  avg=opt.avg,
                  nthreads=opt.nthreads,
                  s_path=ms_input.s_path, t_path=ms_input.t_path)
    command = FdMsInterface.command_from_dict(params)
    if opt.nthreads > 1:
        return command.replace("matching_stats.x", "matching_stats_parallel.x")
    return command


def check_res(res1, res2, avg):
    if avg:
        res1 = float(res1[0])
        res2 = [int(_) for _ in res2[0].split()]
        res2 = sum(res2) / float(len(res2))
        return ["fast (%.4f) != slow (%.4f)" % (res1, res2)
                if round(res1, 4) != round(res2, 4) else False]

    error_lst = []
    for i, (o1, o2) in enumerate(zip(res1, res2)):
        for j, (a, b) in enumerate(zip(o1.split(), o2.split())):
            if a != b:
                error_lst.append(("[%d] %s %s" % (j, a, b)))
    return error_lst


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.verbose else logging.INFO)

    for pref in opt.prefixes:
        ms_input = InputPair.basedir_form(opt.base_dir, pref)
        LG.info("testing on input %s", ms_input)

        res1 = get_output(fast(opt, ms_input))
        res2 = get_output(slow(opt.slow_load_dir, ms_input))

        err_lst = check_res(res1, res2, opt.avg)
        if any(err_lst):
            for err in err_lst:
                if err:
                    print("\t" + err)
        print
        if any(err_lst):
            raise ValueError("ended with errors")


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")

    arg_parser.add_argument("--slow_load_dir", type=str, default=None,
                            help="load results of slow program")
    arg_parser.add_argument("--nthreads", type=int, default=1,
                            help="run parallel version")

    for k in ('double_rank', 'lazy_wl', 'rank_fail', 'use_maxrep_rc', 'use_maxrep_vanilla',
              'load_cst', 'load_maxrep',
              'lca_parents', 'avg', 'answer'):
        args, kwargs = FdMsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    sys.exit(main(arg_parser.parse_args()))
