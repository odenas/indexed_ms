#!/usr/bin/env python


"""
compute matching statistics
"""

import logging
import sys
import argparse

logging.basicConfig(level=logging.INFO)
LG = logging.getLogger(__name__)


def main(opt):
    with open(opt.input) as fd:
        if not opt.no_header:
            fd.readline()
        with open(opt.output, 'w') as ofd:
            for line in fd:
                ofd.write(line.replace('N','').rstrip())


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (denas@adobe.com)")
    arg_parser.add_argument("--input", type=str, default='/dev/stdin',
                            help="input")
    arg_parser.add_argument("--output", type=str, default='/dev/stdout',
                            help="output")
    arg_parser.add_argument("--no_header", action='store_true', default=False,
                            help="Fa file has header")
    sys.exit(main(arg_parser.parse_args()))
