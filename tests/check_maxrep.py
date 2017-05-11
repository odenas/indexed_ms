#!/usr/bin/env python
"""
Run fd_ms on a bunch of inputs while recording space/time usage
"""

import logging
import sys
import argparse

from utils import verbose_args, FullIndex, MaxrepInterface, MsInput, get_output

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def main(opt):
    logging.getLogger().setLevel(logging.DEBUG if opt.v else logging.INFO)

    for i, pref in enumerate(opt.prefixes):
        ispec = MsInput.basedir_form(opt.base_dir, pref)
        command = MaxrepInterface.command_from_dict(dict(s_path=ispec.s_path,
                                                         load_cst=opt.load_cst,
                                                         answer=True))
        res = get_output(command)[0].split()

        with open(ispec.s_path) as fd:
            sidx = FullIndex(fd.read().strip()[::-1])

        query = sidx.string
        checked = set()

        for si in range(len(query)):
            for sj in range(si + 1, len(query)):
                substr = query[si:sj]
                if not sidx.is_node(substr):
                    continue
                i, j = sidx.sa_interval(substr, FullIndex.FWD)
                if (i, j) in checked:
                    continue

                is_max = sidx.is_maximal_s(substr)
                LG.debug("s[%d:%d]=%s, I=[%d, %d], is_max = %s", si, sj, substr, i, j, is_max)
                bin_is_max = (res[i] == res[j - 1] == '1')

                assert is_max == bin_is_max
                checked |= set([(i,j)])


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+",
                            help="input prefixes")

    for k in ('load_cst',):
        args, kwargs = MaxrepInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)
    verbose_args(arg_parser)
    sys.exit(main(arg_parser.parse_args()))
