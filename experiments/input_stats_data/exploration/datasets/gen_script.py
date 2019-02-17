#!/usr/bin/env python


"""
compute matching statistics
"""


import sys
import argparse
from itertools import product
import json

from mstat.dataset import InputPair


def load_json(path):
    with open(path) as fd:
        return json.load(fd)


def create_cmd(opt, config, rep_mut, sim_mut, alp):
    cmd_str = (("python {exec_path} "
                "--rep_blocks {rep_block} "
                "--sim_mut {sim_mut_} "
                "--rep_mut {rep_mut_} "
                "--seed {seed} "
                "--odir {odir_} "
                "rep_sim {alp_} {len_s} {len_t}")
               .format(sim_mut_=sim_mut,
                       rep_mut_=rep_mut,
                       alp_=alp,
                       exec_path=opt.script,
                       odir_=opt.odir,
                       **config))
    return cmd_str


def rename_cmd(opt, config, rep_mut, sim_mut, _alp, which):
    base_path = "rep_{len_s}s_sim_{len_t}t_{_alp}".format(_alp=_alp, **config)
    i = InputPair.basedir_form(opt.odir, base_path)
    opath = getattr(i, which+'_path')
    npath = base_path + "__{rep_mut}_{sim_mut}.{which}".format(**locals())
    return ("mv -v {opath} {npath}".format(**locals()))


def main(opt):
    config = load_json(opt.config)
    param_iter = product(config['rep_mut'],
                         config['sim_mut'],
                         config['alp'])
    for rep_mut, sim_mut, alp in param_iter:
        print(create_cmd(opt, config, rep_mut, sim_mut, alp))
        print(rename_cmd(opt, config, rep_mut, sim_mut, alp, 's'))
        print(rename_cmd(opt, config, rep_mut, sim_mut, alp, 't'))
        print()
    print("touch %s/done" % opt.odir)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Olgert Denas (gertidenas@gmail.com)")
    arg_parser.add_argument('config', type=str, help='config file.')
    arg_parser.add_argument('odir', type=str, help='output dir.')
    arg_parser.add_argument('--script', type=str,
                            default='../../../datasets/generate_input3.py',
                            help='exec')
    sys.exit(main(arg_parser.parse_args()))
