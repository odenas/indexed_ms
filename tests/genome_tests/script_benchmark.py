#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs while recording space/time usage
"""

import logging
import sys
import argparse
import os
import subprocess
import glob

import numpy as np
import pandas as pd

from utils import MsInterface, verbose_args

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def aggregate_max(in_paths_lst, out_path):
    exec_path = ("/Users/denas/Library/Developer/Xcode/DerivedData/"
                 "fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/vec_agg")
    in_paths = " ".join(in_paths_lst)
    command = "{exec_path} {in_paths}".format(**locals())

    with open(out_path, 'w') as ofd:
        LG.info("dumping to %s", out_path)
        return subprocess.check_call("{exec_path} {in_paths}".format(**locals()), shell=True, stdout=ofd)


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
    base_dir = os.path.dirname(opt.t_path)
    out_paths = []
    for i, s_path in enumerate(glob.glob(base_dir + "/" + opt.species + "*.juststring")):
        out_path = "%s__%s" % (opt.t_path.replace(".juststring", ""),
                               os.path.basename(s_path.replace(".juststring", "")))
        command = MsInterface.ms_command_from_dict(dict(s_path=s_path,
                                                        out_path = out_path,
                                                        **base_dict))
        res = subprocess.check_call(command, shell=True)
        out_paths.append(out_path)
    aggregate_max(out_paths, opt.t_path.replace(".juststring", "") + "__" + opt.species)


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
