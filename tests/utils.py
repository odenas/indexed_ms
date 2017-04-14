import logging
import subprocess
import random
import os
from collections import namedtuple, OrderedDict

import numpy as np

LG = logging.getLogger(__name__)

_base_dir_ = ("/Users/denas/Library/Developer/Xcode/DerivedData/"
              "fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/")
PATHS = {"ohleb_fdms" : "%s/fd_ms" % _base_dir_,
         "sada_fdms"  : "%s/sada_fd_ms" % _base_dir_}


def bwt(s):
    from collections import Counter
    s = s + "#"
    A = np.array([list(s[i:] + s[:i]) for i in range(len(s))])
    As = np.vstack(sorted(A, key=tuple))
    CF = Counter(As[:, 0])
    S = sorted(CF)
    C = [0] * (len(CF) + 1)

    for i in range(len(S)):
        if i == 0:
            continue
        C[i] = C[i-1] + CF[S[i-1]]
    C[i+1] = len(s)
    return "".join(As[:, len(s) - 1]), "".join(S), "-".join(map(str, C))


class MsInput(namedtuple('msinput_pair', 's_path, t_path')):
    pass


def _random_input_type((t_path, t_len), (s_path, s_len), source):
    def dump_r(path, l, alp=source):
        with open(path, 'w') as fd:
            i = 0
            while i < l:
                fd.write(random.choice(alp))
                i += 1

    dump_r(t_path, t_len)
    dump_r(s_path, s_len)


def _file_input_type((t_path, t_len), (s_path, s_len), source):
    def _dump_n_chars(ifd, ofd, n):
        i = 0
        while i < n:
            c = ifd.read(1)
            while c == '\n':
                c = infd.read(1)
            ofd.write(c)
            i += 1

    with open(source) as infd:
        with open(t_path, 'w') as outfd:
            _dump_n_chars(infd, outfd, t_len)
        with open(s_path, 'w') as outfd:
            _dump_n_chars(infd, outfd, s_len)


def _mutation_input_type((t_path, t_len), (s_path, s_len), source, mutation_period):
    if t_len > s_len:
        long_path, short_path, = t_path, s_path
        long_len, short_len = t_len, s_len
    else:
        long_path, short_path = s_path, t_path
        long_len, short_len = s_len, t_len

    with open(long_path, 'w') as fd:
        i = 0
        while i < long_len:
            fd.write(random.choice(source))
            i += 1
    import shutil
    shutil.copyfile(long_path, short_path)
    with open(short_path, 'r+') as fd:
        fd.truncate(short_len)
        pos = 0
        while pos < short_len:
            fd.seek(pos)
            fd.write(random.choice(source))
            pos += mutation_period


def create_ms_input(input_type, (t_path, t_len), (s_path, s_len), source, *args):
    def check_len(path):
        res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
        return int(res.split()[0])

    seed = random.randint(0, 10000000)
    LG.warning("SEED: %d", seed)
    random.seed(seed)

    if input_type == "random":
        _random_input_type((t_path, t_len), (s_path, s_len), source)
    elif input_type == "file":
        _file_input_type((t_path, t_len), (s_path, s_len), source)
    elif input_type == "mutation":
        _mutation_input_type((t_path, t_len), (s_path, s_len), source, *args)
    else:
        raise AttributeError("unknown input_type %s" % input_type)

    LG.info("created input of type %s", input_type)
    LG.info("\t%s of length %d", t_path, check_len(t_path))
    LG.info("\t%s of length %d", s_path, check_len(s_path))


class MsInterface(object):
    FDMS_PATH = ("/Users/denas/Library/Developer/Xcode/DerivedData/"
                 "fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/"
                 "fd_ms")
    # name: (required, type, default, help)
    params = OrderedDict(
             s_path = (True, str, None, 'Path of the S string.'),
             t_path = (True, str, None, 'Path of the T string i.e., query'),
             out_path = (False, lambda s: "0" if s is None else str(s), None, 'Dump the ms sdsl::bitvector here.'),
             # sada = (False, bool, False, "Use Sadakane's CST"),
             lazy_wl = (False, bool, False, 'Use lazy weiner links'),
             space_usage = (False, bool, False, 'Report space usage.'),
             time_usage = (False, bool, False, 'Report time usage.'),
             answer = (False, bool, False, 'Print answer.'),
             load_cst = (False, bool, False, 'Load CST of S and Srev.'),
             nthreads = (False, int, 1, 'nr. of threads'),
             runs_progress = (False, int, 0, 'progress msgs for RUNS'),
             ms_progress = (False, int, 0, 'progress msgs for MS'))

    @classmethod
    def as_argparse_kwds(cls, key):
        req, tp, df, h = cls.params[key]
        args = (('' if req else '--') + key,)
        if tp is bool:
            kwargs = {'action': 'store_' + ('false' if df else 'true'),
                      'default': df, 'help': h}
        else:
            kwargs = {'type': tp, 'default': df, 'help': h}
        return (args, kwargs)

    @classmethod
    def as_fdms_params(cls, key, v):
        tp = cls.params[key][1]
        val = (int if tp is bool else tp)(v)
        return "-{key} {val}".format(**locals())


    @classmethod
    def ms_command_from_dict(cls, opt):
        parts = []
        for key, v in opt.iteritems():
            if not (key in cls.params):
                continue
            parts.append(cls.as_fdms_params(key, v))
        # fdms_path = PATHS[('sada_fdms' if opt.get('sada', False) else 'ohleb_fdms')]
        return "{exec_path} {args}".format(exec_path=cls.FDMS_PATH,
                                           args=" ".join(parts))


def verbose_args(arg_parser):
    arg_parser.add_argument("--v", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("--vv", action='store_true',
                            default=False, help="very verbose")
