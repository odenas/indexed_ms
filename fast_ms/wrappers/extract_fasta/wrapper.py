
import logging

import fastaparser


logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


with open(str(snakemake.input)) as fd:
    reader = fastaparser.Reader(fd, parse_method='quick')
    log.info("Parsing %s ...", fd.name)
    for i, seq in enumerate(reader):
        if i != snakemake.params.index:
            continue
        with open(str(snakemake.output), 'wt') as ofd:
            ofd.write(seq.sequence)
