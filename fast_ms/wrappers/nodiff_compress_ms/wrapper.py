
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


tables_path = int(_exists("tables"))


cmd = (
    "{params.exec_path} "
    "-ms_path {snakemake.input.ms}  "
    "-threshold {params.threshold} "
    "-nZeros {params.nzeros} "
    "-nOnes {params.nones} "
    "-negative {params.negative} "
    "-greedy {params.greedy} "
)
log.info(format(cmd))

shell(cmd)

