"""
A Wavelet Tree implementation
"""

from collections import namedtuple
from operator import itemgetter


class Node(namedtuple('node', 'name, parent, children, seq')):
    @property
    def node_symbols(self):
        return sorted(list(set(self.seq)))

    @property
    def mid_symbol(self):
        alp = self.node_symbols
        return alp[len(alp) / 2]

    @property
    def is_leaf(self):
        return len(self.node_symbols) < 2

    @property
    def bin_sequence(self):
        return map(self.char2bin, self.seq)

    def char2bin(self, c):
        return int(c >= self.mid_symbol)

    def bin_rank(self, i, patt):
        s = sum(self.bin_sequence[:i])
        return s if patt == 1 else i - s

    def bin_select(self, cnt, patt):
        assert patt in (0, 1)
        for i, c in enumerate(self.bin_sequence):
            if c == patt:
                cnt -= 1
            if cnt == 0:
                break
        return i

    def bin_select_at_dist(self, i, cnt, patt):
        return self.bin_select(self.bin_rank(i + 1, patt) + cnt, patt)

    def child_seq(self, right):
        return "".join(map(itemgetter(1),
                           filter(lambda t: t[0] == int(right),
                                  zip(self.bin_sequence, self.seq))))


class WTree(object):
    def __init__(self, node_list):
        self.nodes = node_list
        self.root = self.nodes[0]
        self.seq = self.root.seq

    def depth(self, v):
        d = 0
        while v.name != v.parent:
            v = self.nodes[v.parent]
            d += 1
        return d

    def path(self, c):
        p = []
        n = self.root
        while not n.is_leaf:
            i = n.char2bin(c)
            p.append(i)
            n = self.nodes[n.children[i]]
        return p

    def node_of_path(self, p):
        p = list(reversed(p))
        v = self.root
        while p:
            v = self.nodes[v.children[p.pop()]]
        return v

    @classmethod
    def build(cls, pool):
        assert len(pool) == 1
        qool = []
        i = 0
        while pool:
            n = pool.pop()
            if not n.is_leaf:
                l = Node(i + 1, n.name, [], n.child_seq(0))
                r = Node(i + 2, n.name, [], n.child_seq(1))
                n = n._replace(children=(l.name, r.name))
                pool.insert(0, l)
                pool.insert(0, r)
                i += 2
            qool.append(n)
        return cls(qool)

    def rank(self, c, i):
        result = i
        p = list(reversed(self.path(c)))
        v = self.root

        while (not v.is_leaf) and result > 0:
            right = p.pop()
            result = v.bin_rank(result, right)
            v = self.nodes[v.children[right]]
        return result

    def select(self, c, cnt):
        p = self.path(c)
        v = self.node_of_path(p[:-1])
        result = cnt

        while p:
            right = p.pop()
            result = v.bin_select(result, right) + 1
            v = self.nodes[v.parent]
        return result - 1

    def select_at_dist(self, c, i, cnt):
        C = self.seq[i]
        p = self.path(c)
        v = self.root
        i_lst = [i]

        for ck in p[:-1]:
            ik = v.bin_rank(i_lst[-1], ck)
            i_lst.append(ik)
            v = self.nodes[v.children[ck]]

        print i_lst, ik, v, p, C

        if not (C in v.seq):
            i_lst[-1] -= 1
        print [], " ---> [%d]" % 2, v, i_lst[-1], cnt, p[-1]
        jk = v.bin_select_at_dist(i_lst[-1], cnt, p[-1])
        j_lst = [jk]
        v = self.nodes[v.parent]

        for k in reversed(range(len(p) - 1)):
            if not (C in v.seq):
                i_lst[k] -= 1
            print (j_lst, " ---> [%d]" % k, v,
                   i_lst[k], j_lst[-1] - i_lst[k+1], p[k])
            j = v.bin_select_at_dist(i_lst[k], j_lst[-1] - i_lst[k+1], p[k])
            j_lst.append(j)
            v = self.nodes[v.parent]
        result = j_lst[-1]

        exp_result = self.select(c, self.rank(c, i + 1) + cnt)
        assert exp_result == result, "%d != %d" % (exp_result, result)
        return result
