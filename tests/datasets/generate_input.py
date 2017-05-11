#!/usr/bin/env python


"""
generate s and t strings
"""


import logging
import sys
import argparse
import random
import subprocess


LG = logging.getLogger(__name__)


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


def _mutation_input_type((t_path, t_len), (s_path, s_len),
                         source, mutation_period):
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


def create_ms_input(input_type, (t_path, t_len), (s_path, s_len),
                    source, *args):
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


def main(opt):
    create_ms_input(opt.input_type,
                    (opt.t_path, opt.t_len), (opt.s_path, opt.s_len),
                    opt.source, opt.mutation_period)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("input_type",
                            choices=('random', 'file', 'mutation'),
                            help="how to generate.")
    arg_parser.add_argument("source", type=str,
                            help="alphabet for random/mutation, in_file else")

    arg_parser.add_argument('s_path', type=str, help="The S string")
    arg_parser.add_argument("s_len", type=int, help="length of t")
    arg_parser.add_argument('t_path', type=str,
                            help="The T string i.e. the query")
    arg_parser.add_argument("t_len", type=int, help="length of s")

    arg_parser.add_argument("--mutation_period", type=int, default=5,
                            help="mutation period")
    sys.exit(main(arg_parser.parse_args()))
