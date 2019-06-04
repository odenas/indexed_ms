#!/usr/bin/env python


"""
compute matching statistics using a brute force algo
"""


import logging
import random
import sys
import argparse
import subprocess
from typing import Iterator
from tempfile import NamedTemporaryFile as ntp
import numpy as np

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


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


def mstat_py(s, t):
    def iter_prefixes(t: str, i: int) -> Iterator[str]:
        """
        Prefixes of t[i:], from longest to shortest
        """

        for j in reversed(range(i+1, len(t) + 1)):
            yield t[i:j]

    def ms(t: str, s: str, i: int) -> int:
        """
        the longest prefix of t[i:] that occurs in s
        """

        for prefix in iter_prefixes(t, i):
            if prefix in s:
                return len(prefix)
        return 0

    res = []
    for i in range(len(t)):
        res.append(str(ms(t, s, i)))
    return " ".join(res)


def mstat_cpp(s, t):
    exe = "../bin/matching_stats_slow.x"
    with ntp() as s_fd, ntp() as t_fd:
        s_fd.write(s)
        t_fd.write(t)
        s_fd.flush()
        t_fd.flush()

        command = "%s -s_path %s -t_path %s" % (exe, s_fd.name, t_fd.name)
        return get_output(command)[0].strip()


def rnd_str(alp, size, prob=None):
    return "".join(random.choices(population=list(alp), k=size)).encode('ascii')


def main(opt):
    for i in range(opt.repeat):
        s = rnd_str(opt.alp, 10)
        t = rnd_str(opt.alp, 20)

        py = mstat_py(s, t)
        cpp = mstat_cpp(s, t)

        if py != cpp:
            raise ValueError("mismatch: py(%s) != cpp(%s)" % (py, cpp))


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Olgert Denas (gertidenas@gmail.com)"
    )

    arg_parser.add_argument("--repeat", type=int, default=5,
                            help="repeat nr.")
    arg_parser.add_argument("--alp", type=set, default='abc',
                            help="alphabet")
    sys.exit(main(arg_parser.parse_args()))
