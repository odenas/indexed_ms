"""
utilities for dataset generation
"""

import subprocess
import collections
import os
import random
import shutil
import logging

import numpy as np


LG = logging.getLogger("mstat.dataset")


def _check_len(path):
    """
    count the chars in path
    """

    res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
    return int(res.split()[0])


def leave1out(alp):
    """
    create a dict {c: s for c in s} with s = alp - c
    """
    return {c: "".join(set(alp) - set(c)) for c in alp}


def dump_rnd(fd, n, alp, prob=None):
    """
    Dump a random string with the given (char) distribution. Sample usage:

    def t2():
        path = 'dump_rnd.txt'
        np.random.seed(10012)
        alp = list('abcdef')
        prob = np.array([2, 0.1, 0.9, 1, 1, 1], dtype=np.float32) / len(alp)

        with open(path, 'w') as fd:
            dump_rnd(fd, 54, alp, prob)
        print _check_len(path)
        with open(path) as fd:
            cnt = Counter(fd.read().rstrip())
        n = float(sum(cnt.values()))
        return pd.DataFrame([(k, v, v/n) for k, v in cnt.iteritems()], 
                            columns=['c', 'cnt', 'freq'])

    (t2().set_index('c')[['freq']] / 0.15).plot(kind='bar')    
    """

    block_size = 10000
    
    for b in xrange(n / block_size):
        for c in np.random.choice(alp, size=block_size, replace=True, p=prob):
            fd.write(c)

    if n % block_size == 0:
        return

    for c in np.random.choice(alp, size=n % block_size, replace=True, p=prob):
        fd.write(c)


def rnd_textfile(path, text_len, char_counts):
    """
    Create a text file of given length in path. Specify
    symbol frequencies as a counter dict {'char': count}. Sample usage:

    rnd_text('a', 100, 
             Counter('aaa'
                     'bb'
                     'ccc'
                     'dddddd'))
    for k, v in Counter(open('a').read().strip()).iteritems():
        print k, v
    """

    # normalize counts
    alp, prob = zip(*[(k, v / float(sum(char_counts.values())))
                       for k, v in char_counts.iteritems()])

    with open(path, 'w') as fd:
        dump_rnd(fd, text_len, alp, prob)

    LG.info("file %s of length %d", path, _check_len(path))


def mutate_block(path, k, start, end, alp):
    """
    apply k mutations in the block [start, end).

    Sample usage::

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
    mut_positions = sorted(random.sample(xrange(start, end), k))  # w/o replacement

    with open(path, 'r+') as fd:
        for pos in mut_positions:
            fd.seek(pos)
            old_c = fd.read(1)
            new_c = random.choice(i_alp[old_c])

            fd.seek(pos)
            fd.write(new_c)
    return


def mutated_textfile(orig_path, mut_path, num_mutations, alp):
    """
    create a new mutated text file from the given one
    """

    shutil.copyfile(orig_path, mut_path)
    mutate_block(mut_path, num_mutations, 0, _check_len(mut_path), alp)
    LG.info("created %s of length %d from %s of length %d",
            mut_path, _check_len(mut_path), orig_path, _check_len(orig_path))
    return


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