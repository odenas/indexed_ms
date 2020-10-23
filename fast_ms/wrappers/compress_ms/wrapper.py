
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


load_cst = int(_exists('fwd_cst') and _exists('rev_cst'))
load_maxrep = int(_exists("maxrep"))


cmd = (
    "{params.exec_path} "
    "-ms_path {snakemake.input.ms}  "
    "-compression {params.compression}"
)
log.info(format(cmd))

shell(cmd)

