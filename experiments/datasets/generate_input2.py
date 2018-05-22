#!/usr/bin/env python

"""
generate s and t strings
"""


import logging
import os
import sys
import argparse
import random
import subprocess
from collections import namedtuple

LG = logging.getLogger(__name__)


class RepSpec(namedtuple('sr', 'seed_len, seed_nr')):
    def seed_pool(self, alp):
        seed_pool = []
        for seed_idx in range(self.seed_nr):
            curr_seed = [random.choice(alp) for i in range(self.seed_len)]
            seed_pool.append("".join(curr_seed))
        return seed_pool


_input_type_fields = 'dir, alp, slen, stype, tlen, ttype, rep_spec, mperiod'


class InputType(namedtuple('it', _input_type_fields)):
    @property
    def _path(self):
        return os.path.join(self.dir,
                            ("{stype}_{slen}s_"
                             "{ttype}_{tlen}t_"
                             "{alp}_sim{mperiod}")
                            .format(**self._asdict()))

    @property
    def tpath(self):
        return self._path + ".t"

    @property
    def spath(self):
        return self._path + ".s"

    def _dump_rnd(self, path, n):
        with open(path, 'w') as fd:
            i = 0
            while i < n:
                fd.write(random.choice(self.alp))
                i += 1

    def __dump_rep(self, path, n, seed_str_pool):
        with open(path, 'w') as fd:
            char_dumped = 0
            while char_dumped < n:
                seed_str = random.choice(seed_str_pool)
                m = len(seed_str)

                replen = min(random.choice(xrange(1, m)), n - char_dumped)
                repeat_start_position = random.choice(xrange(m - replen))
                assert replen + repeat_start_position <= m

                start_idx = repeat_start_position
                end_idx = min(repeat_start_position + replen, n)
                to_write = seed_str[start_idx:end_idx]

                char_dumped += len(to_write)
                fd.write(to_write)

    def _dump_rep(self, path, n):
        pass

    def _dump_sim(self, from_p, from_l, to_p, to_l, mutation_period):
        import shutil
        shutil.copyfile(from_p, to_p)  # from -> to
        with open(to_p, 'r+') as fd:
            fd.truncate(to_l)
            pos = 0
            while pos < to_l:
                fd.seek(pos)
                fd.write(random.choice(self.alp))
                pos += mutation_period

    def dump(self):
        stypes = ('rnd', 'rep')
        ttypes = ('sim', 'dis')

        assert self.ttype in ttypes
        assert self.stype in stypes

        tp = (self.stype, self.ttype)

        if tp == ('rnd', 'sim'):
            assert self.tlen <= self.slen
            self._dump_rnd(self.spath, self.slen)
            self._dump_sim(self.spath, self.slen, self.tpath, self.tlen,
                           self.mperiod)
        elif tp == ('rnd', 'dis'):
            self._dump_rnd(self.spath, self.slen)
            self._dump_rnd(self.tpath, self.tlen)
        elif tp == ('rep', 'sim'):
            if self.slen > self.tlen:
                self._dump_rep(self.spath, self.slen,
                               self.rep_spec.seed_pool(self.alp))
                self._dump_sim(self.spath, self.slen, self.tpath, self.tlen,
                               self.mperiod)
            else:
                self._dump_rep(self.tpath, self.tlen,
                               self.rep_spec.seed_pool(self.alp))
                self._dump_sim(self.tpath, self.tlen, self.spath, self.slen,
                               self.mperiod)
        else:
            self._dump_rep(self.spath, self.slen,
                           self.rep_spec.seed_pool(self.alp))
            self._dump_rnd(self.tpath, self.tlen)


def main(opt):
    def check_len(path):
        res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
        return int(res.split()[0])

    seed = random.randint(0, 10000000)
    LG.warning("SEED: %d", seed)
    random.seed(seed)

    ispec = InputType(opt.idir, opt.alp,
                      opt.s_len, opt.input_type.split("_")[0],
                      opt.t_len, opt.input_type.split("_")[1],
                      RepSpec(opt.seed_length, opt.seed_pool),
                      opt.mutation_period)
    ispec.dump()

    LG.info("created input of type %s", opt.input_type)
    LG.info("\t%s of length %d", ispec.tpath, check_len(ispec.tpath))
    LG.info("\t%s of length %d", ispec.spath, check_len(ispec.spath))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("input_type",
                            choices=('rnd_sim',
                                     'rnd_dis',
                                     'rep_sim',
                                     'rep_dis'),
                            help="how to generate.")
    arg_parser.add_argument("alp", type=str, help="alphabet")
    arg_parser.add_argument("s_len", type=int, help="len of S")
    arg_parser.add_argument("t_len", type=int, help="len of T, the query")

    arg_parser.add_argument("--idir", type=str, default=".",
                            help="directory to save files.")
    arg_parser.add_argument("--mutation_period", type=int, default=5,
                            help="mutation period")
    arg_parser.add_argument("--seed_length", type=int, default=1000,
                            help="length of seed strings")
    arg_parser.add_argument("--seed_pool", type=int, default=4,
                            help="how many seed strings")
    sys.exit(main(arg_parser.parse_args()))
