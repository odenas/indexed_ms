"""
utilities for dataset generation
"""

import subprocess
import collections
import os
import random
import shutil
import logging
from typing import Iterable, Mapping

import numpy as np

from mstat.interface import get_output


LG = logging.getLogger("mstat.dataset")


def _check_len(path: str) -> int:
    """
    count the chars in path
    """

    res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
    return int(res.split()[0])


def leave1out(alp: Iterable[str]) -> Mapping[str, str]:
    """
    create a dict {c: s for c in s} with s = alp - c
    """
    return {c: "".join(set(alp) - set(c)) for c in alp}


def dump_rnd(fd, n: int, alp: Iterable[str], prob=None) -> None:
    """
    Dump a random string of length `n` from alphabet `alp` with char
    distribution `prob` into file `fd`. Sample usage:

    def t2():
        path = 'dump_rnd.txt'
        np.random.seed(10012)
        alp = list('abcdef')
        prob = np.array([2, 0.1, 0.9, 1, 1, 1], dtype=np.float32) / len(alp)
        assert len(alp) == prob.shape[0]

        with open(path, 'w') as fd:
            dump_rnd(fd, 54, alp, prob)
        print _check_len(path)
        with open(path) as fd:
            cnt = Counter(fd.read().rstrip())
        n = float(sum(cnt.values()))
        return pd.DataFrame([(k, v, v/n) for k, v in cnt.items()],
                            columns=['c', 'cnt', 'freq'])

    (t2().set_index('c')[['freq']] / 0.15).plot(kind='bar')
    """

    block_size = 10000

    for _ in range(n // block_size):
        for c in np.random.choice(alp, size=block_size, replace=True, p=prob):
            fd.write(c)

    if n % block_size == 0:
        return

    for c in np.random.choice(alp, size=n % block_size, replace=True, p=prob):
        fd.write(c)


def rnd_textfile(path: str, text_len: int, char_counts: Mapping[str, int]) -> None:
    """
    Create a text file of given length in path. Specify
    symbol frequencies as a counter dict {'char': count}. Sample usage:

    rnd_textfile('a', 100, Counter('aaa' 'bb' 'ccc' 'dddddd'))
    for k, v in Counter(open('a').read().strip()).items():
        print k, v
    """

    # normalize counts
    alp, prob = zip(*[(k, v / sum(char_counts.values()))
                      for k, v in char_counts.items()])

    with open(path, 'w') as fd:
        dump_rnd(fd, text_len, alp, prob)

    LG.info("file %s of length %d", path, _check_len(path))


def mutate_block(path: str, k: int, start: int, end: int, alp: Iterable[str]) -> None:
    """
    apply k mutations in the block [start, end).

    Sample usage:

    import shutil
    def m1():
        p = 'dump_rnd.txt'
        shutil.copyfile(p, p + '.m')
        mutate_block(p + '.m', 10,  0, 10, list('abcdef'))

        with open(p) as fd:
    orig = fd.read().strip()
    with open(p + '.m') as fd:
        mut = fd.read().strip()

    for i, (a, b) in enumerate(zip(orig, mut)):
        print "%2d" % i,
        if a == b:
            print a
        else:
            print a, "-->", b
    """

    i_alp = leave1out(alp)
    # no replacement
    mut_positions = sorted(random.sample(range(start, end), k))

    with open(path, 'r+') as fd:
        for pos in mut_positions:
            fd.seek(pos)
            old_c = fd.read(1)
            new_c = random.choice(i_alp[old_c])

            fd.seek(pos)
            fd.write(new_c)
    return


def mutated_textfile(orig_path: str, mut_path: str,
                     num_mutations: int, alp: Iterable[str]) -> None:
    """
    create a new mutated text file from the given one
    """

    shutil.copyfile(orig_path, mut_path)
    mutate_block(mut_path, num_mutations, 0, _check_len(mut_path), alp)
    LG.info("created %s of length %d from %s of length %d",
            mut_path, _check_len(mut_path), orig_path, _check_len(orig_path))
    return


def rep_textfile(path: str, n: int, blocks: int, mut: int, alp: Mapping[str, int]):
    """
    generate strings s_1, ..., s_blocks with s_i (i > 1) obtained
    by applying `mut` mutations to s_1.

    dump their concatenation to `path`.
    """

    assert n % blocks == 0
    seed_str_path = path + ".seed"
    mut_str_path = path + ".seed.mutated"
    block_size = n // blocks

    # backup s1 into a file
    rnd_textfile(seed_str_path, block_size, alp)

    # dump s_1
    shutil.copyfile(seed_str_path, path)
    for i in range(0, n, block_size):
        if i == 0:
            continue
        mutated_textfile(seed_str_path, mut_str_path, mut, alp.keys())
        # append s_i
        get_output("cat %s >>%s" % (mut_str_path, path))
        # cleanup
        get_output("rm -f %s" % mut_str_path)
    # remove seed
    get_output("rm -f %s" % seed_str_path)


class InputPair(collections.namedtuple('ip', ['s_path', 't_path'])):
    """
    a pair of files representing a matching statistics input
    """

    @classmethod
    def basedir_form(cls, base_dir, prefix):
        """
        one way to specify inputs
        """
        return cls(os.path.join(base_dir, prefix + ".s"),
                   os.path.join(base_dir, prefix + ".t"))

    @property
    def s_len(self):
        return _check_len(self.s_path)

    @property
    def t_len(self):
        return _check_len(self.t_path)
