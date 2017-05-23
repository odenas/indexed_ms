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


def _repeat_input_type((t_path, t_len), (s_path, s_len), source,
                       seed_len, seed_nr):
    seed_pool = []
    for seed_idx in range(seed_nr):
        current_seed = "".join([random.choice(source) for i in range(seed_len)])
        seed_pool.append(current_seed)

    def dump_str(fd, L, seed_pool=seed_pool):
        char_dumped = 0
        while char_dumped < L:
            seed_str = seed_pool[random.choice(xrange(len(seed_pool)))]
            a, b = random.sample(xrange(len(seed_str)), 2)
            l, h = min(a, b), max(a, b)
            to_write = (seed_str[l:h])[0:(L - char_dumped)]
            char_dumped += len(to_write)
            fd.write(to_write)

    with open(t_path, 'w') as fd:
        dump_str(fd, t_len)
    with open(s_path, 'w') as fd:
        dump_str(fd, s_len)


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


INPUT_TYPES = {'rnd': (_random_input_type, []),
               'file': (_file_input_type, []),
               'mut': (_mutation_input_type, ['mutation_period']),
               'rep': (_repeat_input_type, ['seed_length', 'seed_pool'])
               }


def main(opt):
    def check_len(path):
        res = subprocess.check_output("/usr/bin/wc -c %s" % path, shell=True)
        return int(res.split()[0])

    seed = random.randint(0, 10000000)
    LG.warning("SEED: %d", seed)
    random.seed(seed)

    input_function, arg_names = INPUT_TYPES[opt.input_type]
    input_function((opt.t_path, opt.t_len), (opt.s_path, opt.s_len), opt.source, *[getattr(opt, n) for n in arg_names])

    LG.info("created input of type %s", input_type)
    LG.info("\t%s of length %d", t_path, check_len(t_path))
    LG.info("\t%s of length %d", s_path, check_len(s_path))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")

    arg_parser.add_argument("input_type", choices=INPUT_TYPES,
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
    arg_parser.add_argument("--seed_length", type=int, default=1000,
                            help="mutation period")
    arg_parser.add_argument("--seed_pool", type=int, default=4,
                            help="mutation period")
    sys.exit(main(arg_parser.parse_args()))
