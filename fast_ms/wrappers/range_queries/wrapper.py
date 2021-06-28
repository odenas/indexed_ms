
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


ridx_opts = ""
if _exists("ridx"):
    ridx_opts = format("-ridx_path {snakemake.input.ridx}")


cmd = (
    "{params.exec_path} "
    "-ms_path {snakemake.input.ms} "
    "{ridx_opts} "
    "-from_idx {params.from_idx} "
    "-to_idx {params.to_idx} "
    "-block_size {params.block_size} "
    "-compression {params.compression} "
    "-algo {params.algo} "
    "-op {params.op} "
    " >{snakemake.output}"
)
log.info(format(cmd))

shell(cmd)

