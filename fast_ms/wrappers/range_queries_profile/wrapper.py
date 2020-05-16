
import os
from pathlib import Path
from snakemake.shell import shell
from snakemake.utils import format
import logging



params = snakemake.params
exe = Path(str(params.bin_dir), "range_queries_profile.x")

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(Path(exe).stem)

def _exists(which):
    try:
        return Path(str(getattr(snakemake.input, which))).is_file()
    except AttributeError:
        return False


ridx_opts = ""
if _exists("ridx"):
    ridx_opts = format("-ridx_path {snakemake.input.ridx}")


cmd = (
    "{exe} "
    "-ms_path {snakemake.input.ms} "
    "-compression {params.compression} "
    "-block_size {params.block_size} {ridx_opts} "
    "-algo {params.algo} -op {params.op} "
    "-from_max_idx {params.from_max_idx} -range_size {params.range_size} -niter {params.niter} "
    " >{snakemake.output}"
)
log.info(format(cmd))

shell(cmd)

