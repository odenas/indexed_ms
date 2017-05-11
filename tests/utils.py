import logging
import subprocess
import os
from collections import namedtuple, OrderedDict

import numpy as np

LG = logging.getLogger(__name__)

_base_dir_ = ("/Users/denas/Library/Developer/Xcode/DerivedData/"
              "fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/")
PATHS = {"ohleb_fdms": "%s/fd_ms" % _base_dir_,
         "sada_fdms": "%s/sada_fd_ms" % _base_dir_}


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
    @classmethod
    def basedir_form(cls, base_dir, prefix):
        return cls(os.path.join(base_dir, prefix + ".s"),
                   os.path.join(base_dir, prefix + ".t"))


class XcodeBinaryInterface(object):
    @classmethod
    def _check_class_members(cls):
        for name in ("params", "EXEC_PATH"):
            if not hasattr(cls, name):
                msg = "Class should define a `%s` classmember" % name
                raise AttributeError(msg)

    @classmethod
    def strfarg(cls, key, val):
        """
        Stringify the key, val pair to an argument of the executable. Fetch the
        type of the argument from cls.params
        """

        cls._check_class_members()

        arg_type = cls.params[key][1]
        if arg_type is bool:  # bools are encoded as ints
            arg_type = int
        return "-{key} {val}".format(val=arg_type(val), key=key)

    @classmethod
    def command_from_dict(cls, opt):
        """
        Generate a command that can be run on a terminal
        from the given dictionary `opt`.
        """

        cls._check_class_members()

        parts = []
        for key, v in opt.iteritems():
            if not (key in cls.params):
                continue
            parts.append(cls.strfarg(key, v))
        return "{exec_path} {args}".format(exec_path=cls.EXEC_PATH,
                                           args=" ".join(parts))

    @classmethod
    def as_argparse_kwds(cls, key):
        """generate args and kwargs to pass to argparse.add_arguments()"""

        cls._check_class_members()

        req, tp, df, h = cls.params[key]
        args = (('' if req else '--') + key,)
        if tp is bool:
            kwargs = {'action': 'store_' + ('false' if df else 'true'),
                      'default': df, 'help': h}
        else:
            kwargs = {'type': tp, 'default': df, 'help': h}
        return (args, kwargs)


class MsInterface(XcodeBinaryInterface):
    EXEC_PATH = os.path.join(_base_dir_, "fd_ms")

    # name: (required, type, default, help)
    params = OrderedDict(
             s_path=(True, str, None, 'Path of the S string.'),
             t_path=(True, str, None, 'Path of the T string i.e., query'),
             out_path=(False, lambda s: "0" if s is None else str(s), None,
                       'Dump the ms sdsl::bitvector here.'),
             rank_fail=(False, bool, False, "Use the rank-and-fail strategy."),
             lazy_wl=(False, bool, False, 'Use lazy Weiner links'),
             use_maxrep=(False, bool, False, 'maxrep vector for Weiner links'),
             space_usage=(False, bool, False, 'Report space usage.'),
             time_usage=(False, bool, False, 'Report time usage.'),
             answer=(False, bool, False, 'Print answer.'),
             load_cst=(False, bool, False, 'Load CST of S and Srev.'),
             nthreads=(False, int, 1, 'nr. of threads'),
             runs_progress=(False, int, 0, 'progress msgs for RUNS'),
             ms_progress=(False, int, 0, 'progress msgs for MS'))


def verbose_args(arg_parser):
    arg_parser.add_argument("--v", action='store_true',
                            default=False, help="verbose")
    arg_parser.add_argument("--vv", action='store_true',
                            default=False, help="very verbose")


def get_output(command):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")
