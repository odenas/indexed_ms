"""
Commandline interfaces to the C++ binaries
"""

import logging
import sys
import subprocess
import os
import argparse
from collections import OrderedDict
from typing import Mapping, Tuple, Any, Callable

LG = logging.getLogger(__name__)
ArgumentSpecType = Tuple[bool, Callable, Any, str]

_base_dir_ = '/home/odenas/code/matching_statistics/indexed_ms/fast_ms/bin'
sys.stderr.write("*** using base directory: %s ***\n" % _base_dir_)


def get_output(command: str) -> str:
    """
    run and parse output

    :param str command:
    :return: list of output lines
    """

    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, shell=True).decode()
    LG.debug("got: " + res)
    return res.strip().split("\n")


def default_arg_parser(doc, add_verbose=True, add_seed=True):
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
        arg_parser.add_argument("-v", "--verbose",
                                action='store_true',
                                default=False,
                                help="verbose")
    if add_seed:
        arg_parser.add_argument("--seed", type=int, default=None,
                                help="Seed of random number generator.")

    return arg_parser


class CommandLineArguments(object):
    """
    Utility to translate between binary arguments and python's argparse.

    This is meant to be subclassed -- subclass must define members
     - params
     - EXEC_PATH
    """

    params: Mapping[str, ArgumentSpecType]
    EXEC_PATH: str

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
        Generate a command that calls a c++ binary from
        the given dictionary `opt`.

        :param dict opt: dictionary key -> val arguments
        """

        cls._check_class_members()

        parts = []
        for key, val in opt.items():
            if key not in cls.params:
                continue
            parts.append(cls.strfarg(key, val))
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
    """
    interface to the dump_cst binary
    """

    EXEC_PATH = os.path.join(_base_dir_, "dump_cst.x")

    # name: (required, type, default, help)
    params = OrderedDict([CommonArgumentSpecs.s_path])


class FdMsInterface(CommandLineArguments):
    """
    Interface to the fd_ms binary
    """

    EXEC_PATH = os.path.join(_base_dir_, "matching_stats.x")
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
        ('avg', (False, bool, False, 'Print average MS.')),
        CommonArgumentSpecs.load_cst,
        ('nthreads', (False, int, 1, 'nr. of threads')),
        ('nslices', (False, int, 1, 'nr. of blocks'))
    ])
