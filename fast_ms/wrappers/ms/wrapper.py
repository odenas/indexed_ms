
import os
from pathlib import Path
from snakemake.shell import shell
from snakemake.utils import format
import logging



params = snakemake.params

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(Path(params.exec_path).stem)

def _exists(which):
    try:
        return Path(str(getattr(snakemake.input, which))).is_file()
    except AttributeError:
        return False

def _try_remove(file_path):
    file_path = Path(file_path)

    if file_path.is_file():
        file_path.unlink()
        if not file_path.is_file():
            log.info("* removed %s", file_path)


load_cst = int(_exists('fwd_cst') and _exists('rev_cst'))
load_maxrep = int(_exists("maxrep"))
if len(snakemake.output) != 1:
    raise AttributeError("Expecting a single output, the .ms file, got %d instead: %s",
                         len(snakemake.output), str(snakemake.output))


cmd = (
    "{params.exec_path} "
    "-s_path {snakemake.input.s} -t_path {snakemake.input.t} "
    "-lca_parents {params.lca_parents} -rank_fail {params.rank_fail} -double_rank {params.double_rank} -lazy_wl {params.lazy_wl} "
    "-load_cst {load_cst} -load_maxrep {load_maxrep} "
    "-nthreads {snakemake.threads} "
)
log.info(format(cmd))

shell(cmd)

if hasattr(snakemake.params, "remove_temp") and snakemake.params.remove_temp:
    file_id = str(snakemake.output)[:-3]
    log.info("removing temporary files starting with %s ...", file_id)
    _try_remove(file_id + ".runs")

    for i in range(snakemake.threads):
        _try_remove(file_id + (".ms.%d" % i))
        _try_remove(file_id + (".runs.%d" % i))
        _try_remove(file_id + (".runs.%d_%d" % (i, i + 1)))

if hasattr(snakemake.params, "remove_idx") and snakemake.params.remove_idx:
    file_id = str(snakemake.input.s)
    _try_remove(file_id + ".fwd.stree")
    _try_remove(file_id + ".rev.stree")
    _try_remove(file_id + ".rev.maxrep")

