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

    @staticmethod
    def dump(fn, s):
        LG.debug("\tdumping to %s ..." % fn)
        with open(fn, 'w') as fd:
            fd.write(s)
        return fn

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

    @staticmethod
    def get_bp(in_fname, out_fname, bp_exec="./bp"):
        command = "{bp_exec} {in_fname} > {out_fname}".format(**locals())
        res = subprocess.check_output(command, shell=True)
        return res.strip()

    @classmethod
    def random_build(cls, opt):
        prefix = "random"
        def getp(i):
            fname = str(prefix) + cls.ft2fn(cls._fields[i])
            return os.path.join(opt.base_dir, fname)

        l, L = opt.l, max(opt.l, opt.L)
        for prefix in range(opt.n):
            t = "".join([random.choice(opt.alphabet)
                         for i in range(random.randint(l, L))])
            s = "".join([random.choice(opt.alphabet)
                         for i in range(random.randint(l, L))])

            LG.info("creating input %s ..." % prefix)
            c0 = cls.dump(getp(0), t)
            c1 = cls.dump(getp(1), s)
            c2 = cls.get_bp(getp(1), getp(2))
            s = s[::-1]
            c3 = cls.dump(getp(3), s)
            c4 = cls.get_bp(getp(3), getp(4))
            return cls._make([c0, c1, c2, c3, c4])

    @classmethod
    def file_build(cls, in_file, len_t, len_s, prefix, base_dir):
        def getp(i):
            fname = prefix + cls.ft2fn(cls._fields[i])
            return os.path.join(base_dir, fname)

        with open(in_file) as fd:
            cnt = 0
            txt = []
            for line in fd:
                for c in line.rstrip():
                    if cnt >= len_t + len_s:
                        break
                    txt.append(c)
                    cnt += 1
        assert len(txt) == len_t + len_s

        def dump(fn, i, j, reverse):
            with open(fn, 'w') as fd:
                if reverse:
                    k = j - 1
                    while k >= i:
                        fd.write(txt[k])
                        k -= 1
                else:
                    k = i
                    while k < j:
                        fd.write(txt[k])
                        k += 1
                fd.write("\n")
            return fn

        c0 = dump(getp(0), 0, len_t, False)

        c1 = dump(getp(1), len_t, len(txt), False)
        c2 = cls.get_bp(getp(1), getp(2))

        c3 = dump(getp(3), len_t, len(txt), True)
        c4 = cls.get_bp(getp(3), getp(4))

        return cls._make([c0, c1, c2, c3, c4])


def main(opt):
    if opt.base_text == 'random':
        for i in range(opt.n):
            Input.random_build(opt)
    else:
        prefix = os.path.basename(opt.in_file) + ("%d_%d" % (opt.len_t, opt.len_s))
        Input.file_build(opt.in_file, opt.len_t, opt.len_s, prefix, opt.base_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("base_text", choices=('random', 'non_random'),
                            help="a filename or 'random'")
    arg_parser.add_argument("--base_dir", type=str, default="./input_data",
                            help="base dir")

    nrgrp = arg_parser.add_argument_group('options non-random data')
    nrgrp.add_argument("--in_file", type=str, help="base text")
    nrgrp.add_argument("--len_t", type=int, default=10,
                       help="t is the prefix of the input")
    nrgrp.add_argument("--len_s", type=int, default=5,
                       help="s follows t in the input")

    rgrp = arg_parser.add_argument_group('options random data')
    rgrp.add_argument("--n", type=int, help="how many string pairs")
    rgrp.add_argument("--l", type=int, default=7, help="length lbound")
    rgrp.add_argument("--L", type=int, default=10, help="length hbound")
    rgrp.add_argument("--alphabet", type=str, default='ab', help="alphabet")
    sys.exit(main(arg_parser.parse_args()))
