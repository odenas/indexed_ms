#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import subprocess
import sys
import random
import argparse
from collections import Counter

import numpy as np


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def bwt(s):
    s = s + "#"
    A = np.array([list(s[i:] + s[:i]) for i in range(len(s))])
    As = np.vstack(sorted(A, key = tuple))
    CF=Counter(As[:,0])
    S = sorted(CF)
    C = [0] * (len(CF) + 1)

    for i in range(len(S)):
        if i == 0:
            continue
        C[i] = C[i-1] + CF[S[i-1]]
    C[i+1] = len(s)
    return "".join(As[:,len(s) - 1]), "".join(S), "-".join(map(str, C))


def get_string(n, alp):
    alp = list(alp)
    return "".join([random.choice(alp) for i in range(n)])

def get_bp(s, bp_exec = "./bp"):
    command = "echo {s} | {bp_exec} /dev/stdin".format(**locals())
    res = subprocess.check_output(command, shell=True)
    return res.strip()


def main(opt):
    for pair in range(opt.n):
        len_t = random.randint(opt.l, max(opt.l, opt.L))
        len_s = random.randint(opt.l, max(opt.l, opt.L))

        t = get_string(len_t, ['a', 'b'])
        s = get_string(len_t, set(t))
        while set(t) != set(s):
            s = get_string(len_t, set(t))
        assert set(s) == set(t)

        s_bp = get_bp(s)
        srev = s[::-1]
        srev_bp = get_bp(srev)
        print t, s, s_bp, srev, srev_bp


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("n", type=int, help="how many string pairs")
    arg_parser.add_argument("--l", type=int, default=7, help="string length lower bound")
    arg_parser.add_argument("--L", type=int, default=10, help="string length higher bound")
    sys.exit(main(arg_parser.parse_args()))
