import logging
import subprocess
import random
import os
from collections import namedtuple


LG = logging.getLogger(__name__)


def bwt(s):
    import numpy as np
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


II = namedtuple('ii', ('t_path, s_fwd_path, s_fwd_bp_path, '
                       's_rev_path, s_rev_bp_path'))


class Input(II):
    alphabet = ['a', 'b']

    @classmethod
    def ft2fn(cls, ft):
        return ft.replace('_path', '.txt')

    @classmethod
    def fn2ft(cls, ft):
        return ft.replace('.txt', '_path')

    @staticmethod
    def dump(fn, s):
        LG.debug("\tdumping to %s ..." % fn)
        with open(fn, 'w') as fd:
            fd.write(s)
        return fn

    @classmethod
    def _parse(cls, prefix, base_dir=".", check_exist=True):
        def make_p(ft):
            path = os.path.join(base_dir, prefix + cls.ft2fn(ft))
            if check_exist and (not os.path.exists(path)):
                raise ValueError("bad path. %s does not exist" % (path,))
            LG.debug("\tparsed %s", path)
            return path

        LG.debug("parsing input %s ..." % prefix)
        return cls._make(map(make_p, cls._fields))

    @staticmethod
    def get_bp(in_fname, out_fname, bp_exec="./bp"):
        command = "{bp_exec} {in_fname} > {out_fname}".format(**locals())
        res = subprocess.check_output(command, shell=True)
        return res.strip()

    @classmethod
    def random_build(cls, opt):
        def getp(i):
            fname = str(prefix) + cls.ft2fn(cls._fields[i])
            return os.path.join(opt.base_dir, fname)

        l, L = opt.l, max(opt.l, opt.L)
        for prefix in range(opt.n):
            t = "".join([random.choice(opt.alphabet)
                         for i in range(random.randint(l, L))])
            s = "".join([random.choice(opt.alphabet)
                         for i in range(random.randint(l, L))])

            LG.info("creating input %s ..." % prefix)
            c0 = cls.dump(getp(0), t)
            c1 = cls.dump(getp(1), s)
            c2 = cls.get_bp(getp(1), getp(2))
            s = s[::-1]
            c3 = cls.dump(getp(3), s)
            c4 = cls.get_bp(getp(3), getp(4))
        return cls._make([c0, c1, c2, c3, c4])

    @classmethod
    def file_build(cls, in_file, len_t, len_s, prefix, base_dir):
        def getp(i):
            fname = prefix + cls.ft2fn(cls._fields[i])
            return os.path.join(base_dir, fname)

        with open(in_file) as fd:
            cnt = 0
            txt = []
            for line in fd:
                for c in line.rstrip():
                    if cnt >= len_t + len_s:
                        break
                    txt.append(c)
                    cnt += 1
        assert len(txt) == len_t + len_s

        def dump(fn, i, j, reverse):
            with open(fn, 'w') as fd:
                if reverse:
                    k = j - 1
                    while k >= i:
                        fd.write(txt[k])
                        k -= 1
                else:
                    k = i
                    while k < j:
                        fd.write(txt[k])
                        k += 1
                fd.write("\n")
            return fn

        c0 = dump(getp(0), 0, len_t, False)

        c1 = dump(getp(1), len_t, len(txt), False)
        c2 = cls.get_bp(getp(1), getp(2))

        c3 = dump(getp(3), len_t, len(txt), True)
        c4 = cls.get_bp(getp(3), getp(4))

        return cls._make([c0, c1, c2, c3, c4])


class MsCommand(object):
    @classmethod
    def fast(self, input_dir, prefix,
             space_usage, mem_usage,
             time_usage,
             answer, verb,
             path_to_exec):
        cmd_templ = ("{exec_path} -d {dir} -p {prefix} "
                     "-S {mach_mem} -s {sp} -t {tm} -a {ans} -v {verb}")
        return (cmd_templ
                .format(exec_path=path_to_exec,
                        dir=input_dir,
                        prefix=prefix,
                        sp=int(space_usage),
                        mach_mem = int(mem_usage),
                        tm=int(time_usage),
                        ans=int(answer),
                        verb=int(verb)))

    @classmethod
    def slow(self, input_dir, prefix, path_to_exec):
        return ("python {exec_path} {dir} {prefix}"
                .format(exec_path=path_to_exec,
                        dir=input_dir, prefix=prefix))


def get_output(command, *args, **kwargs):
    LG.debug("running: " + str(command))
    res = subprocess.check_output(command, *args, shell=True, **kwargs)
    LG.debug("got: " + res)
    return res.strip().split("\n")
