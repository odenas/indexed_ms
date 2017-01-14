import logging
import subprocess
import random
import os
from collections import namedtuple

import numpy as np

LG = logging.getLogger(__name__)
FDMS_PATH = ("/Users/denas/Library/Developer/Xcode/DerivedData/"
             "fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms")


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


def _random_input_type(t_path, len_t, s_path, len_s, alphabet):
    def dump_r(path, l, alp=alphabet):
        with open(path, 'w') as fd:
            i = 0
            while i < l:
                fd.write(random.choice(alp))
                i += 1

    dump_r(t_path, len_t)
    dump_r(s_path, len_s)


def _file_input_type(t_path, len_t, s_path, len_s, infile):
    def _dump_n_chars(ifd, ofd, n):
        i = 0
        while i < n:
            c = ifd.read(1)
            while c == '\n':
                c = infd.read(1)
            ofd.write(c)
            i += 1

    with open(infile) as infd:
        with open(t_path, 'w') as outfd:
            _dump_n_chars(infd, outfd, len_t)
        with open(s_path, 'w') as outfd:
            _dump_n_chars(infd, outfd, len_s)


def _mutation_input_type(t_path, len_t, s_path, len_s, alphabet, mt_period):
    if len_t > len_s:
        long_path, short_path, long_len, shorlen_t = (t_path, s_path, len_t, len_s)
    else:
        long_path, short_path, long_len, shorlen_t = (s_path, t_path, len_s, len_t)

    with open(long_path, 'w') as fd:
        i = 0
        while i < long_len:
            fd.write(random.choice(alphabet))
            i += 1
    import shutil
    shutil.copyfile(long_path, short_path)
    with open(short_path, 'r+') as fd:
        fd.truncate(shorlen_t)
        pos = 0
        while pos < shorlen_t:
            fd.seek(pos)
            fd.write(random.choice(alphabet))
            pos += mt_period


def create_input(input_spec, len_t, len_s, source, input_type,
                 mutation_period):
    seed = random.randint(0, 10000000)
    LG.warning("SEED: %d", seed)
    random.seed(seed)
    if input_type == "random":
        _random_input_type(input_spec.t_path, len_t, input_spec.s_path, len_s, source)
    elif input_type == "file":
        _file_input_type(input_spec.t_path, len_t, input_spec.s_path, len_s, source)
    elif input_type == "mutation":
        _mutation_input_type(input_spec.t_path, len_t, input_spec.s_path, len_s, source, mutation_period)
    else:
        raise AttributeError("unknown input type %s" % input_type)

    LG.info("created input wrt %s", str(input_spec))
    LG.info("%s of length %d", input_spec.t_path, input_spec.len_t)
    LG.info("%s of length %d", input_spec.s_path, input_spec.len_s)


class InputSpec(namedtuple('iii', 'base_dir, prefix')):
    def _path(self, which):
        return os.path.join(self.base_dir, "%s%s.txt" % (self.prefix, which))

    @property
    def t_path(self): return self._path("t")

    @property
    def s_path(self): return self._path("s")


    @property
    def paths(self):
        return [self.t_path, self.s_path]

    @property
    def len_t(self):
        res = subprocess.check_output("/usr/bin/wc -c %s" % self.t_path, shell=True)
        return int(res.split()[0])

    @property
    def len_s(self):
        res = subprocess.check_output("/usr/bin/wc -c %s" % self.s_path, shell=True)
        return int(res.split()[0])


class MsCommand(object):
    @classmethod
    def fast(self, input_spec,
             lazy_wl, sada_st,
             space_usage, time_usage,
             answer,
             verb,
             runs_prg, ms_prg,
             path_to_exec):
        cmd_templ = ("{exec_path} "
                     "-d {dir} -p {prefix} "
                     "-l {lazy} -sada {sada} "
                     "-s {sp} -t {tm} "
                     "-a {ans} "
                     "-v {verb} "
                     "-runs_progress {runs_prg} "
                     "-ms_progress {ms_prg} ")
        return (cmd_templ
                .format(exec_path=path_to_exec,
                        dir=input_spec.base_dir, prefix=input_spec.prefix,
                        lazy=int(lazy_wl), sada=int(sada_st),
                        sp=int(space_usage), tm=int(time_usage),
                        ans=int(answer),
                        verb=int(verb),
                        runs_prg=int(runs_prg), ms_prg=int(ms_prg),
                        ))

    @classmethod
    def slow(self, input_spec, path_to_exec):
        return ("python {exec_path} {dir} {prefix}"
                .format(exec_path=path_to_exec,
                        dir=input_spec.base_dir, prefix=input_spec.prefix))


def get_output(command, *args, **kwargs):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, *args, shell=True, **kwargs)
    LG.debug("got: " + res)
    return res.strip().split("\n")
