#!/usr/bin/env python

"""
Check correctness of maxrep
"""

import logging
import sys

from bin_interfaces import DumpMaxrepInterface, MsInput, get_output, default_arg_parser
from indexes import FullIndex

#from utils import verbose_args, FullIndex, MaxrepInterface, MsInput, get_output

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger()


def main(opt):
    """
    main thred of program
    """

    logging.getLogger().setLevel(logging.DEBUG if opt.verbose else logging.INFO)

    for i, pref in enumerate(opt.prefixes):
        ispec = MsInput.basedir_form(opt.base_dir, pref)
        command = DumpMaxrepInterface.command_from_dict(dict(s_path=ispec.s_path,
                                                             load_cst=opt.load_cst,
                                                             txt_format=True))
        res = get_output(command)[0].split()

        with open(ispec.s_path) as fd:
            sidx = FullIndex(fd.read().strip()[::-1])

        for ((si, sj), substr, (i, j), is_max) in sidx.maxrep_iter():
            LG.debug("s[%d:%d]=%s, I=[%d, %d], is_max = %s",
                     si, sj, substr, i, j, is_max)
            bin_is_max = (res[i] == res[j - 1] == '1')

            assert is_max == bin_is_max


if __name__ == "__main__":
    arg_parser = default_arg_parser(__doc__)

    arg_parser.add_argument("base_dir", type=str, help="base dir")
    arg_parser.add_argument("prefixes", type=str, nargs="+", help="input prefixes")

    for k in ('load_cst',):
        args, kwargs = DumpMaxrepInterface.as_argparse_kwds(k)
        arg_parser.add_argument(*args, **kwargs)

    sys.exit(main(arg_parser.parse_args()))
