from itertools import product
from pathlib import PosixPath as Path
import re

# wrt to the Snakefile
bin_dir = Path("../../bin").absolute()
wrapper_dir = Path("../../wrappers").absolute()
assert bin_dir.is_dir()
idir = "../input"
assert Path(idir).is_dir()


# defines input generation in the ./input directory

def strpiid(input_id, __iid_re=r"([a-z]{3})_(\d+)s_([a-z]{3})_(\d+)t_([a-z]+)"):
    m = re.match(__iid_re, input_id)
    if m is None:
        _msg = f"Couldn't parse {input_id} with pattern {__iid_re}"
        raise ValueError(_msg)
    raw_tup = m.groups()
    return raw_tup[0], int(raw_tup[1]), raw_tup[2], int(raw_tup[3])


def strfiid(it, il, qt, ql, alp):
    return f"{it}_{il}s_{qt}_{ql}t_{alp}"


input_index_types = ["rnd", "rep"]
input_index_lengths = (300,)
input_query_types = ["dis", "sim"]
input_query_lengths = (50, 150)
alphabet = "abcde"
_iid_iter = product(input_index_types, input_index_lengths,
                    input_query_types, input_query_lengths, [alphabet])
input_ids = [strfiid(*tup) for tup in _iid_iter]


ms_par = Path(bin_dir, "matching_stats_parallel.x")
ms_ = Path(bin_dir, "matching_stats.x")
ms_slow = Path(bin_dir, "matching_stats_slow.x")
split_ms = Path(bin_dir, "split_ms.x")
print_int_ms = Path(bin_dir, "print_int_ms.x")
print_bin_ms = Path(bin_dir, "print_ms.x")
print_freq = Path(bin_dir, "print_freq.x")
dump_maxrep = Path(bin_dir, "dump_maxrep.x")
dump_cst = Path(bin_dir, "dump_cst.x")
range_query = Path(bin_dir, "range_queries.x")
range_index = Path(bin_dir, "dump_range_index.x")
compress_ms = Path(bin_dir, "compress_ms.x")

for p in [ms_par, ms_slow, split_ms,
          print_int_ms, dump_maxrep, dump_cst]:
    assert p.is_file(), "{pp} not found".format(pp=p)

iids = [str(s.stem) for s in Path(idir).glob("*.s")]


class ipair():
    avail_compr = ('none', 'rrr', 'delta',
                   'succint', 'nibble', 'rle',
                   )

    def __init__(self, pair_id):
        self.pair_id = pair_id

        self.s = self.pair_id + ".s"
        self.t = self.pair_id + ".t"

        self.rev_cst = self.s + ".rev.stree"
        self.fwd_cst = self.s + ".fwd.stree"

        self.maxrep = self.s + ".rev.maxrep"

        self.mstat = self.pair_id + ".mstat"
        self.ms_path = self.s + "_" + self.t + ".ms"
        self.freq_path = self.s + "_" + self.t + ".freq"

    def ridx(self, block_size):
        return self.ms_path + ".none." + str(block_size) + ".ridx"

    def compr_ms_path(self, compr):
        if compr == "{compr}":
            return self.ms_path + "." + compr

        if compr not in self.avail_compr:
            raise ValueError("bad compression %s. must be one of %s" % (compr, self.avail_compr))
        return self.ms_path + "." + compr
