import os
from pathlib import PosixPath as Path

# wrt to the Snakefile
bin_dir = Path("../../bin").absolute()
assert bin_dir.is_dir()
idir = "../input"
assert Path(idir).is_dir()


ms_par = Path(bin_dir, "matching_stats_parallel.x")
ms_slow = Path(bin_dir, "matching_stats_slow.x")
split_ms = Path(bin_dir, "split_ms.x")
print_int_ms = Path(bin_dir, "print_int_ms.x")
print_bin_ms = Path(bin_dir, "print_ms.x")
dump_maxrep = Path(bin_dir, "dump_maxrep.x")
dump_cst = Path(bin_dir, "dump_cst.x")
range_query = Path(bin_dir, "range_queries.x")
range_index = Path(bin_dir, "dump_range_index.x")
compress_ms = Path(bin_dir, "compress_ms.x")

for p in [ms_par, ms_slow, split_ms,
          print_int_ms, dump_maxrep, dump_cst]:
    assert p.is_file(), "{pp} not found".format(pp=p)

iids   = [str(s.stem) for s in Path(idir).glob("*.s")]

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
        self.ms_path = self.t + ".ms"

    def ridx(self, block_size):
        return self.ms_path + ".none." + str(block_size) + ".ridx"

    def compr_ms_path(self, compr):
        if compr == "{compr}":
            return self.ms_path + "." + compr

        if not compr in self.avail_compr:
            raise ValueError("bad compression %s. must be one of %s" % (compr, self.avail_compr))
        return self.ms_path + "." + compr


