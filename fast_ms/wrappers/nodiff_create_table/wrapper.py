
from pathlib import Path
from snakemake.shell import shell
from snakemake.utils import format
import logging


params = snakemake.params

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(Path(params.exec_path).stem)


out_path = Path(str(snakemake.output))
out_dir = out_path.parent

cmd = (
    "{params.exec_path} "
    "{params.from_threshold} {params.to_threshold} "
    "{params.nzeroes} {params.nones} "
    "{params.negative} "
    "0 "
    "{out_dir} "
)
log.info(format(cmd))
shell(cmd)

if not Path(str(snakemake.output)).is_file():
    # rename
    cmd = (
        "mv "
        "table-nodiff-{params.from_threshold}-{params.nzeroes}-{params.nones}-{params.negative} "
        "{snakemake.output}"
    )
    shell(cmd)
