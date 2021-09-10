
from pathlib import Path
from snakemake.shell import shell
from snakemake.utils import format
import logging


params = snakemake.params

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(Path(params.exec_path).stem)

outdir = Path(str(snakemake.output)).parent

cmd = (
    "{params.exec_path} "
    "{params.from_threshold} {params.to_threshold} "
    "{params.nzeros} {params.nones} "
    "{params.negative} "
    "0 "
    "{outdir} "
)
log.info(format(cmd))
shell(cmd)

shell("touch {snakemake.output}")
