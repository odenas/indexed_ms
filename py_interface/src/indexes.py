
from collections import Counter, namedtuple, OrderedDict
from typing import Iterator, Tuple

import numpy as np
import pandas as pd

# todo: see https://github.com/BenLangmead/comp-genomics-class


class Bwt(namedtuple('bwt', 'text, bwt, alp, cumm_sym_counts, sa')):
    """
    the BWT of a string
    """

    @classmethod
    def from_text(cls, s: str):
        """
        instantiate Bwt from the given text

        :param str s: the text
        :return: Bwt instance
        """

        s += '#'
        a = np.array([list(s[i:] + s[:i]) for i in range(len(s))])
        sa = np.vstack(sorted(a, key=tuple))
        sym_counts = Counter(sa[:, 0])
        alphabet = sorted(sym_counts)
        cummulative_sym_counts = [0] * (len(sym_counts) + 1)

        for i in range(len(alphabet)):
            if i == 0:
                continue
            cummulative_sym_counts[i] = cummulative_sym_counts[i-1] + sym_counts[alphabet[i-1]]
        cummulative_sym_counts[i+1] = len(s)

        return cls(s,
                   "".join(sa[:, len(s) - 1]),
                   "".join(alphabet),
                   "-".join(map(str, cummulative_sym_counts)),
                   sa)

    @property
    def index_table(self) -> str:
        """
        A dataframe with varous indexes of text.
        """

        s, bwt_s, *_, sa = self

        def i_to_sa_suff(i):
            return "".join(sa[i][:sa[i].tolist().index("#") + 1])

        a = OrderedDict([
            ('i', range(len(s))),
            ('s_i', list(s)),
            ('BWT', list(bwt_s)),
            ('SA', [len(s) - sa[i].tolist().index("#") for i in range(sa.shape[0])]),
            ('suff_SA', [i_to_sa_suff(i) for i in range(sa.shape[0])])
        ])
        return pd.DataFrame(a).set_index('i')


class FullIndex(object):
    """
    Index data structures for a given text
    """

    FWD, REV = 0, 1

    def __init__(self, s):
        self.string,  self.bwt_s, self.alp, self.C, self.sa = Bwt.from_text(s)

        self.tabs = {self.FWD: Bwt.from_text(s).index_table,
                     self.REV: Bwt.from_text(s[::-1]).index_table}

    def is_node(self, substr: str) -> bool:
        """
        is `substr` the label of a proper node?
        """

        def _next_char(sa_str):
            l = len(substr)
            h = min(l + 1, len(sa_str))
            assert h <= l + 1
            return sa_str[l:h]

        i, j = self.sa_interval(substr, FullIndex.FWD)
        ext_chars = set(map(_next_char,
                            self.tabs[FullIndex.FWD][i:j].suff_SA))
        return len(ext_chars) > 1

    def count(self, s: str, dir: int) -> int:
        tab = self.tabs[dir]
        pref_bool_idx = tab.suff_SA.apply(lambda suff: suff.startswith(s))
        return pref_bool_idx.sum()

    def sa_interval(self, s: str, dir: int) -> Tuple[int, int]:
        tab = self.tabs[dir]
        pref_bool_idx = tab.suff_SA.apply(lambda suff: suff.startswith(s))
        idx = tab.loc[pref_bool_idx].index.tolist()
        return min(idx), max(idx) + 1

    def is_maximal_s(self, s) -> bool:
        i, j = self.sa_interval(s, self.FWD)
        maxleft = len(set(self.tabs[self.FWD][i:j].BWT)) > 1
        return maxleft and self.is_node(s)

    def maxrep_iter(self) -> Iterator:
        checked = set()

        yield ((0, 0),
               "",
               tuple(self.sa_interval('', FullIndex.FWD)),
               self.is_maximal_s(""))

        for si in range(len(self.string)):
            for sj in range(si + 1, len(self.string)):
                substr = self.string[si:sj]
                if not self.is_node(substr):
                    continue
                i, j = self.sa_interval(substr, FullIndex.FWD)
                if (i, j) in checked:
                    continue
                yield ((si, sj), substr, (i, j), self.is_maximal_s(substr))
                checked |= {(1, 2)}


class MatchingStatistics(namedtuple('ms', ['index', 'query', 'matching_statistics', 'ms', 'runs'])):
    @staticmethod
    def _ms(t: str, s: str, i: int) -> int:
        """the longest prefix of t[i:] that occurs in s"""

        def iter_prefixes(t: str, i: int) -> Iterator[str]:
            "from longest to shortest"
            for j in reversed(range(i+1, len(t) + 1)):
                yield t[i:j]

        for prefix in iter_prefixes(t, i):
            if prefix in s:
                return len(prefix)
        return 0

    @classmethod
    def ms_table(cls, t: str, s: str):
        """
        A dataframe with matching statistics and the intermediate vectors runs, ms
        """

        a = pd.DataFrame(OrderedDict([
            ('i', [-1] + list(range(len(t)))),
            ('t_i', [''] + list(t)),
            ('MS', [1] + list(map(lambda i: cls._ms(t, s, i), range(len(t)))))
        ]))
        a['nzeros'] = [0] + (a[1:].MS.values - a[0:len(t)].MS.values + 1).tolist()
        a['ms'] = [''] + list(map(lambda i: "".join(['0'] * i) + '1', a.nzeros[1:]))
        a['runs'] = [-1] + list(map(lambda s: int(s == '1'), a.ms[1:]))
        a['t_MS'] = [t[_sl[0]:_sl[0]+_sl[1]] for _sl in zip(a.i, a.MS)]
        s_idx = FullIndex(s)
        a['F'] = [s_idx.count(_p, FullIndex.FWD) for _p in a.t_MS]
        f_step = [False] + [a.F[i] != a.F[i-1] for i in range(1, len(a.F))]
        # for frequency we have a modified runs
        runs_freq = a.runs.tolist()
        for i in range(1, a.shape[0]):
            c1 = a.runs.iloc[i]
            c2 = a.F.iloc[i] != a.F.iloc[i - 1]
            if c1 and c2:
                runs_freq[i] = 0
        a['runs_freq'] = runs_freq

        return a.set_index('i')

    @classmethod
    def from_text(cls, index: str, query: str):
        """
        instantiate MatchingStatistics from given input
        :param str index:
        :param str query:
        :return: MatchingStatistics instance
        """

        a = cls.ms_table(query, index)
        return cls(index, query, a.MS.tolist(), a.ms.values.tolist(), a.runs.tolist())


