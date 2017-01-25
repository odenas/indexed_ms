#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs while recording space/time usage
"""

import logging
import sys
import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

from utils import MsInterface, verbose_args

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def get_output(command):
    LG.debug("running: " + str(command))
    with tempfile.NamedTemporaryFile(delete=False) as ofd:
        LG.debug("dumping to : " + str(ofd.name))
        res = subprocess.check_call(command, shell=True, stdout=ofd)
        ofd.flush()
        X = np.loadtxt(ofd.name, dtype=np.uint32)
    LG.debug("got array of size %s ", str(X.shape))
    return X


def get_species_paths(species, base_dir="./genome_data"):
    import glob
    #LG.error("clipping species chromosomes to 3")
    return glob.glob(base_dir + "/" + species + "*.juststring")

def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    base_dict = dict(t_path = opt.t_path,
                     lazy_wl=opt.lazy_wl,
                     sada=False,
                     space_usage=False, time_usage=False,
                     answer=True,
                     load_cst=opt.load_cst,
                     runs_progress=opt.runs_progress,
                     ms_progress=opt.ms_progress)

    npy_paths = []
    for i, s_path in enumerate(get_species_paths(opt.species)):
        command = MsInterface.ms_command_from_dict(dict(s_path=s_path,
                                                        **base_dict))
        npy_paths.append("%s__%s.npy" % (os.path.basename(opt.t_path), os.path.basename(s_path)))
        np.save(npy_paths[-1], get_output(command))
    LG.info("repacking ...")
    X = pd.DataFrame(dict(zip(map(lambda s: s.split("__")[1].replace(".juststring.npy", ""), npy_paths),
                              map(np.load, npy_paths))))
    out_path = "%s__%s.mk.pkl" % (opt.t_path.replace(".juststring", ""), opt.species)
    LG.info("dumping dataframe (%s) to %s", str(X.shape), out_path)
    X.to_pickle(out_path)

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("t_path", type=str, help="The T string input path")
    arg_parser.add_argument("species", choices=("Mus_musculus",
                                                "Homo_sapiens",
                                                "Danio_rerio"),
                            help="The species from where to draw S strings")

    for k in ('lazy_wl', 'runs_progress', 'ms_progress', 'load_cst'):
        args, kwargs = MsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)

    arg_parser.add_argument("--label", type=str, default='default',
                            help="ad a label col to output")
    arg_parser.add_argument("--repeat", type=int, default=1,
                            help="repeat each experiment")
    verbose_args(arg_parser)
    sys.exit(main(arg_parser.parse_args()))
