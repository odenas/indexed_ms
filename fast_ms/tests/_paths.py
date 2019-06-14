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
dump_maxrep = Path(bin_dir, "dump_maxrep.x")
dump_cst = Path(bin_dir, "dump_cst.x")

for p in [ms_par, ms_slow, split_ms,
          print_int_ms, dump_maxrep, dump_cst]:
    assert p.is_file()

iids   = [str(s.stem) for s in Path(idir).glob("*.s")]

class ipair():
    def __init__(self, pair_id):
        self.pair_id = pair_id

        self.s = self.pair_id + ".s"
        self.t = self.pair_id + ".t"

        self.rev_cst = self.s + ".rev.stree"
        self.fwd_cst = self.s + ".fwd.stree"

        self.maxrep = self.s + ".rev.maxrep"

        self.mstat = self.pair_id + ".mstat"
        self.ms_path = self.t + ".ms"
