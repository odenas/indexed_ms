#!/usr/bin/env python


"""
compute matching statistics

/Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms
"""


import logging
import sys
import argparse

from generate_strings import Input

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def bwt(s):
    import numpy as np
    from collections import Counter
    s = s + "#"
    A = np.array([list(s[i:] + s[:i]) for i in range(len(s))])
    As = np.vstack(sorted(A, key=tuple))
    CF = Counter(As[:, 0])
    S = sorted(CF)
    C = [0] * (len(CF) + 1)

    for i in range(len(S)):
        if i == 0:
            continue
        C[i] = C[i-1] + CF[S[i-1]]
    C[i+1] = len(s)
    return "".join(As[:, len(s) - 1]), "".join(S), "-".join(map(str, C))


def iter_prefixes(t, i):
    "from longest to shortest"

    for j in reversed(range(i+1, len(t) + 1)):
        yield t[i:j]


def ms(t, s, i):
    for prefix in iter_prefixes(t, i):
        if prefix in s:
            return len(prefix)
    return 0


def main(opt):
    t = open(getattr(opt, Input._fields[0])).read().strip()
    s = open(getattr(opt, Input._fields[1])).read().strip()
    print getattr(opt, Input._fields[0]), "".join([str(ms(t, s, i)) for i in range(len(t))])


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument(Input._fields[0], type=str, help="t")
    arg_parser.add_argument(Input._fields[1], type=str, help="s")
    arg_parser.add_argument(Input._fields[2], type=str, help="bp of s")
    arg_parser.add_argument(Input._fields[3], type=str, help="reversed s")
    arg_parser.add_argument(Input._fields[4], type=str, help="bp reversed s")
    sys.exit(main(arg_parser.parse_args()))
