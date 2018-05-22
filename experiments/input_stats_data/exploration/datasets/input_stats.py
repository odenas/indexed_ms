#!/usr/bin/env python


"""
compute input statistics
"""


from collections import OrderedDict, namedtuple
from contextlib import contextmanager
import logging
import os
import random
import subprocess
import sys

from mstat.interface import (CommandLineArguments,
                             _base_dir_,
                             default_arg_parser,
                             get_output,
                             CommonArgumentSpecs)
sys.path.append("../../../datasets")
from generate_input3 import leave1out


logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


class InputStatsInterface(CommandLineArguments):
    EXEC_PATH = os.path.join(_base_dir_, "input_stats.x")
    params = OrderedDict([CommonArgumentSpecs.s_path,
                          CommonArgumentSpecs.t_path,
                          CommonArgumentSpecs.load_cst,
                          CommonArgumentSpecs.load_maxrep])


class PathParts(namedtuple('pp', 'alp, st, sl, tt, tl, rep_mut, sim_mut')):
    @classmethod
    def parse_path(cls, s_path):
        a, b = os.path.basename(s_path).split('__')
        st, sl, tt, tl, alp = a.split('_')
        rep_mut, sim_mut = b[:-2].split('_')
        return cls._make((alp, st,
                          int(sl[:-1]), tt, int(tl[:-1]),
                          int(rep_mut), int(sim_mut)))

    @property
    def mut_positions(self):
        return sorted(random.sample(xrange(self.tl), self.sim_mut))


def check_len(path):
    res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
    return int(res.split()[0])


@contextmanager
def get_stats(pparts, s_path, cmd_dict):
    t_path = s_path + "trial.t"
    if os.path.exists(t_path):
        raise AttributeError("t_path (%s) exists" % t_path)

    i_alp = leave1out(pparts.alp)

    import shutil
    shutil.copyfile(s_path, t_path)

    with open(t_path, 'r+') as fd:
        for pos in pparts.mut_positions:
            fd.seek(pos)
            old_c = fd.read(1)
            fd.seek(pos)
            fd.write(random.choice(i_alp[old_c]))
        fd.truncate(pparts.tl)
        LG.info("%s --> %s (length %d, %d mutations)",
                s_path, t_path, check_len(t_path), pparts.sim_mut)

    cmd_dict = dict(cmd_dict, **{'t_path': t_path})
    yield get_output(InputStatsInterface.command_from_dict(cmd_dict))

    LG.info("removing %s", t_path)
    os.remove(t_path)


def dump_res(fd, res, trial, s_path):
    for i, r in enumerate(res):
        if i == 0 and trial == 0:
            fd.write(r + ",trial,b_path\n")
        else:
            fd.write(r + "," +
                     str(trial + 1) + "," +
                     os.path.basename(s_path) + "\n")


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.verbose
                                 else logging.INFO)

    if opt.output != '/dev/stdout' and os.path.exists(opt.output):
        LG.error("output file (%s) exsts. Exiting ...", opt.output)
        return 1

    path_parts = PathParts.parse_path(opt.s_path)
    cmd_dict = dict({'s_path': opt.s_path,
                     'load_cst': opt.load_cst,
                     'load_maxrep': opt.load_maxrep})
    with open(opt.output, 'w') as ofd:
        for i in range(opt.repeat):
            with get_stats(path_parts, opt.s_path, cmd_dict) as res:
                dump_res(ofd, res, i, opt.s_path)


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)

    for k in InputStatsInterface.params:
        args, kwargs = InputStatsInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    arg_parser.add_argument("--repeat", type=int, default=1,
                            help="repeat experiment")
    sys.exit(main(arg_parser.parse_args()))
