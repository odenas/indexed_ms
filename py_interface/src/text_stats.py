
import array
from collections import Counter
import itertools
import os
import random

import numpy as np
import pandas as pd

from mstat.dataset import _check_len


def input_char_table(path):
    """
    a count/frequency table of symbols in the given file
    """

    with open(path) as fd:
        cnt = Counter(itertools.chain.from_iterable(fd))

    n = float(sum(cnt.values()))
    return pd.DataFrame([(k, v, v / n) for k, v in cnt.iteritems()],
                        columns=['c', 'cnt', 'freq'])


def mutation_table(path1, path2, name1=None, name2=None):
    """
    a table of mutated chars in path1 to obtain path2
    """

    if name1 is None:
        name1 = os.path.basename(path1)
    if name2 is None:
        name2 = os.path.basename(path2)

    with open(path1) as fd1:
        with open(path2) as fd2:
            txt1 = itertools.chain.from_iterable(fd1)
            txt2 = itertools.chain.from_iterable(fd2)
            data = [(a, b, a == b) for (a, b) in zip(txt1, txt2)]
    return pd.DataFrame(data, columns=[name1, name2, 'is_mut'])


def infer_alphabet(path):
    """
    get alphabet (as a counter) from given file
    """

    return Counter(dict([(row.c, row.cnt)
                         for i, row in input_char_table(path).iterrows()]))


def repfile_bitmap(path, blocks):
    """
    represent the repeats-file as a bitmap of blocks rows with
    1-values at mutations. bitmap[i,j] == 1 iff the i-th block
    has a mutation (wrt 0-th block) at position j.
    """

    _x = array.array('B')
    with open(path) as fd:
        _x.fromfile(fd, _check_len(path))
    x = np.array(_x, dtype=np.uint8).reshape((blocks, -1))
    return ((x - np.copy(x[0])) != 0)


def block_sample_iter(txt, block_size, sample_size):
    """
    sample random blocks of text from the given text string
    """

    for i in range(sample_size):
        start_idx = random.randint(0, len(txt) - block_size)
        yield start_idx, txt[start_idx:(start_idx + block_size)]
