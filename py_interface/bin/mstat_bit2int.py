"""
Print the MS values by reading them off a (text) ms vector
"""
import argparse
import sys
from typing import Iterator


def bit2int_iter(bit_iter: Iterator[int]) -> Iterator[int]:
    """
    Iterate over MS vaues given a (bit) ms iterator.
    :param iterable[int] bit_iter: Iterator over ms values
    :return: iterable[int]
    """

    prev_ms = -1
    nones = 0
    nzeros = 0
    for i, val in enumerate(bit_iter):
        if val == 1:
            ms_val = prev_ms - nzeros + 1
            assert ms_val == i - (2 * nones)
            yield (ms_val)
            nones += 1
            prev_ms = ms_val
        elif val == 0:
            nzeros += 1
        else:
            raise ValueError("expecting a bit-vec in input. %d-th value = %d" % (i, val))


def main(opt):
    with open(opt.bit_ms) as fd:
        for ms_val in bit2int_iter(map(int, fd)):
            print(ms_val)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog="Olgert Denas (gertidenas@gmail.com)")
    arg_parser.add_argument('bit_ms', type=str, help='Path of the bit ms vector')
    sys.exit(main(arg_parser.parse_args()))
