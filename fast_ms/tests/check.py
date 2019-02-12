#!/usr/bin/env python

"""
Test whether programs produce the same output
"""

import logging
import sys
from difflib import ndiff

from mstat.interface import default_arg_parser, get_output, FdMsInterface
from mstat.dataset import InputPair


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def fast(ms_input, nthreads):
    params = dict(lazy_wl=False, rank_fail=False, double_rank=False, lca_parents=False,
                  use_maxrep_vanilla=False, use_maxrep_rc=False,
                  answer=True, avg=False,
                  nthreads=nthreads, s_path=ms_input.s_path, t_path=ms_input.t_path)
    command = FdMsInterface.command_from_dict(params)
    return command.replace("matching_stats.x", "matching_stats_parallel.x")


def check_line(l1, l2, name):
    print("checking line %s" % name)
    err = [(" " if (v1 == v2) else "*") for (v1, v2) in zip(l1, l2)]
    err = [e + (" " if (len(v2) > 1) else "") for (e, v2) in zip(err, l2)]

    print(" ".join(l1))
    print(" ".join(err))
    print(" ".join(l2))
    print()
    return bool("".join(err).strip())


def check_res(lines1, lines2):
    assert len(lines1) == len(lines2)

    e = []
    for l1, l2 in zip(lines1, lines2):
        e.append(check_line(l1.strip().split(), l2.strip().split(), ""))
    return e


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.verbose else logging.INFO)

    ms_input = InputPair(opt.s_path, opt.t_path)
    LG.info("testing on input %s", ms_input)

    res1 = get_output(fast(ms_input, 1))
    res2 = get_output(fast(ms_input, opt.nthreads))

    if any(check_res(res1, res2)):
        raise ValueError("ended with errors")


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)
    arg_parser.add_argument("s_path", type=str, help="base dir")
    arg_parser.add_argument("t_path", type=str, help="base dir")
    arg_parser.add_argument("nthreads", type=int, help="base dir")
    sys.exit(main(arg_parser.parse_args()))
