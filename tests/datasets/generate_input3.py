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

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")


def leave1out(alp):
    """
    create a dict of the type {c: s} with s = alp - c
    """
    return {c: "".join(set(alp) - set(c)) for c in alp}


class InputType(namedtuple('it', ['dir', 'alp',
                                  'slen', 'stype', 'tlen', 'ttype'])):
    @property
    def _path(self):
        return os.path.join(self.dir,
                            ("{stype}_{slen}s_{ttype}_{tlen}t_{alp}"
                             .format(**self._asdict())))

    @property
    def tpath(self):
        return self._path + ".t"

    @property
    def spath(self):
        return self._path + ".s"

    def _dump_rnd(self, path, n):
        with open(path, 'w') as fd:
            for i in xrange(n):
                fd.write(random.choice(self.alp))

    def _dump_rep(self, path, n, mut, blocks):
        """
        generate `blocks` string s1, ..., sk with si (i > 1) obtained from s1
        by applying `mut` mutations.

        dump their concatenation to `path`.
        """

        assert n % blocks == 0
        seed_str_path = path + ".seed"
        mut_str_path = path + ".seed.mutated"
        block_size = n / blocks
        # backup s1 into a file
        self._dump_rnd(seed_str_path, block_size)

        # dump s1
        get_output("cat %s >>%s" % (seed_str_path, path))
        for i in range(0, n, block_size):
            if i == 0:
                continue
            self._dump_sim(seed_str_path, block_size,
                           mut_str_path, block_size,
                           mut)
            # concatenate si
            get_output("cat %s >>%s" % (mut_str_path, path))

        # cleanup
        get_output("rm -f %s" % (' '.join([mut_str_path, seed_str_path]),))


    def _dump_sim(self, from_p, from_l, to_p, to_l, mutations):
        """
        Generate a string S by applying `mutations` mutations to text T[0:`to_l`].
        Obtain T from `from_p` and dump S[0:`to_l`]
        """

        i_alp = leave1out(self.alp)

        import shutil
        shutil.copyfile(from_p, to_p)  # from -> to
        with open(to_p, 'r+') as fd:
            mut_positions = sorted(random.sample(xrange(to_l), mutations))
            for pos in mut_positions:
                fd.seek(pos)
                old_c = fd.read(1)
                fd.seek(pos)
                fd.write(random.choice(i_alp[old_c]))
            fd.truncate(to_l)
            # LG.info("%s --> %s (%d mutations)", from_p, to_p, mutations)

    def dump(self, sim_mut, rep_mut, rep_blocks):
        assert self.ttype in ('sim', 'dis')
        assert self.stype in ('rnd', 'rep')
        tp = (self.stype, self.ttype)

        if tp == ('rnd', 'sim'):
            assert self.tlen <= self.slen
            self._dump_rnd(self.spath, self.slen)
            self._dump_sim(self.spath, self.slen,
                           self.tpath, self.tlen,
                           sim_mut)
        elif tp == ('rnd', 'dis'):
            self._dump_rnd(self.spath, self.slen)
            self._dump_rnd(self.tpath, self.tlen)
        elif tp == ('rep', 'sim'):
            self._dump_rep(self.spath, self.slen, rep_mut, rep_blocks)
            self._dump_sim(self.spath, self.slen, self.tpath, self.tlen,
                           sim_mut)
        elif tp == ('rep', 'dis'):
            self._dump_rep(self.spath, self.slen, rep_mut, rep_blocks)
            self._dump_rnd(self.tpath, self.tlen)


def main(opt):
    def check_len(path):
        res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
        return int(res.split()[0])

    seed = (random.randint(0, 10000000) if opt.seed is None else opt.seed)
    LG.warning("SEED: %d", seed)
    random.seed(seed)

    ispec = InputType(opt.odir, opt.alp,
                      opt.sizes[0], opt.input_type.split("_")[0],
                      opt.sizes[1], opt.input_type.split("_")[1])
    ispec.dump(opt.sim_mut, opt.rep_mut, opt.rep_blocks)

    LG.info("created input of type %s", opt.input_type)
    LG.info("\t%s of length %d", ispec.tpath, check_len(ispec.tpath))
    LG.info("\t%s of length %d", ispec.spath, check_len(ispec.spath))


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (gertidenas@gmail.com)")

    arg_parser.add_argument("input_type",
                            choices=('rnd_sim',
                                     'rnd_dis',
                                     'rep_sim',
                                     'rep_dis'),
                            help="Type of dataset.")
    arg_parser.add_argument("alp", type=str, help="alphabet")
    arg_parser.add_argument("sizes", type=int, nargs=2,
                            help="index and query lengths")

    arg_parser.add_argument("--odir", type=str, default=".",
                            help="directory to save files.")
    arg_parser.add_argument("--sim_mut", type=int, default=5,
                            help="nr. of mutations for similar string")
    arg_parser.add_argument("--rep_mut", type=int, default=5,
                            help="nr. of mutations for repeat string")
    arg_parser.add_argument("--rep_blocks", type=int, default=5,
                            help="nr. of blocks for repeat string")
    arg_parser.add_argument("--seed", type=int, default=None,
                            help="Seed of random number generator.")
    sys.exit(main(arg_parser.parse_args()))
