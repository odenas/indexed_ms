#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
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


def get_string(n):
    return "".join([random.choice(['a', 'b', 'c']) for i in range(n)])


def main(opt):
    for pair in range(opt.n):
        len_t = random.randint(15, 20)
        len_s = random.randint(15, 20)

        t, s = get_string(len_t), get_string(len_s)
        assert set(s) == set(t)
        bwt_fw, C, S = bwt(s)
        bwt_rev, _, _ = bwt(s[::-1])
        print t, s, bwt_fw, bwt_rev, C, S


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("n", type=int, help="how many string pairs")
    sys.exit(main(arg_parser.parse_args()))
