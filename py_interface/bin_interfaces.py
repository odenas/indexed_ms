"""
Interfaces to the C++ binaries
"""


import logging
import os
import sys
import subprocess
import argparse
from collections import OrderedDict, namedtuple


LG = logging.getLogger(__name__)

_base_dir_candidates_ = [
        "/home/brt/Documents/projects/matching_statistics/indexed_ms/fast_ms/bin",
        ]
_existing_bdirs_ = filter(os.path.exists, _base_dir_candidates_)
_base_dir_ = _existing_bdirs_[0]
print >>sys.stderr, "*** using base directory: %s ***" % _base_dir_


def default_arg_parser(doc, add_verbose=True):
    """
    instantiate an argparser
    :param str doc:
    :param add_verbose:
    :return:
    """
    arg_parser = argparse.ArgumentParser(
            description=doc,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (gertidenas@gmail.com)")
    if add_verbose:
        arg_parser.add_argument("-v", "--verbose", action='store_true', default = False, help="verbose")
    return arg_parser


class MsInput(namedtuple('msinput_pair', 's_path, t_path')):
    @classmethod
    def basedir_form(cls, base_dir, prefix):
        return cls(os.path.join(base_dir, prefix + ".s"),
                   os.path.join(base_dir, prefix + ".t"))


class CommandLineArguments(object):
    """
    Utility to translate between binary arguments and python's argparse.

    This is meant to be subclassed -- subclass must define members
     - params
     - EXEC_PATH
    """

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
        type of the argument from `cls.params`

        :param str key: the name of the argyument
        :param val: the value of the argument
        """

        cls._check_class_members()

        arg_type = cls.params[key][1]
        if arg_type is bool:  # bools are encoded as ints
            arg_type = int
        return "-{key} {val}".format(val=arg_type(val), key=key)

    @classmethod
    def command_from_dict(cls, opt):
        """
        Generate a command that calls a c++ binary from the given dictionary `opt`.

        :param dict opt: dictionary key -> val arguments
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
        """
        generate args and kwargs to pass to argparse.add_arguments()

        :param str key: argument name
        """

        cls._check_class_members()

        req, tp, df, h = cls.params[key]
        args = (('' if req else '--') + key,)
        if tp is bool:
            kwargs = {'action': 'store_' + ('false' if df else 'true'),
                      'default': df, 'help': h}
        else:
            kwargs = {'type': tp, 'default': df, 'help': h}
        return args, kwargs


class CommonArgumentSpecs(object):
    """
    name says it
    """

    # name: (required, type, default, help)
    s_path = ('s_path', (True, str, None, 'Path of the index string.'))
    t_path = ('t_path', (True, str, None, 'Path of the query string.'))
    out_path = ('out_path', (False, lambda s: "0" if s is None else str(s), None, 'Dump results here.'))
    load_cst = ('load_cst', (False, bool, False, 'Load CST of S and reverse(S).'))
    load_maxrep = ('load_maxrep', (False, bool, False, 'Load the maxrep vector.'))


class FdMsInterface(CommandLineArguments):
    """
    Interface to the fd_ms binary
    """

    EXEC_PATH = os.path.join(_base_dir_, "matching_stats.opt.x")
    # name: (required, type, default, help)
    params = OrderedDict([
        CommonArgumentSpecs.s_path,
        CommonArgumentSpecs.t_path,
        CommonArgumentSpecs.out_path,
        ('double_rank', (False, bool, False, "Use double instead of single rank.")),
        ('rank_fail', (False, bool, False, "Use the rank-and-fail strategy.")),
        ('lazy_wl', (False, bool, False, 'Use lazy Weiner links')),
        ('use_maxrep_rc', (False, bool, False, 'maxrep vector for Weiner links with rank&check')),
        ('use_maxrep_vanilla', (False, bool, False, 'maxrep vector for Weiner links')),
        ('lca_parents', (False, bool, False, 'Use lca instead of consecutive parent sequence')),
        CommonArgumentSpecs.load_maxrep,
        ('time_usage', (False, bool, False, 'Report time usage.')),
        ('answer', (False, bool, False, 'Print answer.')),
        CommonArgumentSpecs.load_cst,
        ('nthreads', (False, int, 1, 'nr. of threads'))
    ])


class DumpMaxrepInterface(CommandLineArguments):
    """
    interface to the dump_maxrep binary
    """

    EXEC_PATH = os.path.join(_base_dir_, "dump_maxrep.x")

    # name: (required, type, default, help)
    params = OrderedDict([CommonArgumentSpecs.s_path,
                          ('txt_format', (False, bool, False, "Dump as txt.")),
                          CommonArgumentSpecs.load_cst])



class DumpCstInterface(CommandLineArguments):
    pass


def get_output(command):
    """
    run and parse output

    :param str command:
    :return: list of output lines
    """

    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True)
    LG.debug("got: " + res)
    return res.strip().split("\n")
