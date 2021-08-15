
from pathlib import Path
from snakemake.shell import shell
from snakemake.utils import format
import logging


params = snakemake.params

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(Path(params.exec_path).stem)

outdir = Path(str(snakemake.output)).parent

threshold_minus = int(params.threshold) - 1

cmd = (
    "{params.exec_path} "
    "{params.threshold} "
    "0 "
    "{threshold_minus} "
    "{params.max_n} "
    "{params.max_n} "
    "{params.negative} "
    "0 "
    "{outdir} "
)
log.info(format(cmd))
shell(cmd)

shell("touch {snakemake.output}")
