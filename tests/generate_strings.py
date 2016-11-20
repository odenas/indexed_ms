#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import subprocess
import sys
import random
import argparse
import os
from collections import namedtuple


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


II = namedtuple('ii', ('t_path, s_fwd_path, s_fwd_bp_path, '
                       's_rev_path, s_rev_bp_path'))


class Input(II):
    alphabet = ['a', 'b']

    @classmethod
    def ft2fn(cls, ft):
        return ft.replace('_path', '.txt')

    @classmethod
    def fn2ft(cls, ft):
        return ft.replace('.txt', '_path')

    def read(self):
        dt = []
        for ft in self:
            dt.append(open(getattr(self, ft)).read().strip())
        return dt

    @classmethod
    def _parse(cls, prefix, base_dir=".", check_exist=True):
        def make_p(ft):
            path = os.path.join(base_dir, prefix + cls.ft2fn(ft))
            if check_exist and (not os.path.exists(path)):
                raise ValueError("bad path. %s does not exist" % (path,))
            LG.debug("\tparsed %s", path)
            return path

        LG.debug("parsing input %s ..." % prefix)
        return cls._make(map(make_p, cls._fields))

    @classmethod
    def _build(cls, len_t, len_s, prefix, base_dir="."):
        def get_bp(s, bp_exec="./bp"):
            command = "echo {s} | {bp_exec} /dev/stdin".format(**locals())
            res = subprocess.check_output(command, shell=True)
            return res.strip()

        def dump(fn, data_iter):
            LG.debug("\tdumping to %s ..." % fn)
            with open(fn, 'w') as fd:
                for c in data_iter:
                    fd.write(c)
                fd.write("\n")
            return fn

        def getp(i):
            fname = prefix + cls.ft2fn(cls._fields[i])
            return os.path.join(base_dir, fname)

        LG.info("creating input %s ..." % prefix)
        c0 = dump(getp(0),
                  (random.choice(cls.alphabet) for i in range(len_t)))
        s = "".join([random.choice(cls.alphabet) for i in range(len_s)])
        c1 = dump(getp(1), s)
        c2 = dump(getp(2), get_bp(s))
        s = s[::-1]
        c3 = dump(getp(3), s)
        c4 = dump(getp(4), get_bp(s))

        return cls._make([c0, c1, c2, c3, c4])


def main(opt):
    for i in range(opt.n):
        len_t = random.randint(opt.l, max(opt.l, opt.L))
        len_s = random.randint(opt.l, max(opt.l, opt.L))
        Input._build(len_t, len_s, str(i), opt.base_dir)

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("n", type=int, help="how many string pairs")
    arg_parser.add_argument("--l", type=int, default=7, help="length lbound")
    arg_parser.add_argument("--L", type=int, default=10, help="length hbound")
    arg_parser.add_argument("--base_dir", type=str, default="input_data", help="base dir")
    sys.exit(main(arg_parser.parse_args()))
