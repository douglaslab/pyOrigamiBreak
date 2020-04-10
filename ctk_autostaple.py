#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-05-21 14:35:43
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import argparse
import sys


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i",   "--input",    type=str,
                        help="Cadnano json file")

    parser.add_argument("-rule",   "--rule",     type=str, default='',
                        help="Staple rule")

    parser.add_argument("-seq",   "--sequence", type=str, default=None,
                        help="Sequence file in txt")

    parser.add_argument("-nsol",   "--nsol",     type=int,
                        help="Number of solutions", default=10)

    parser.add_argument("-v",   "--verbose",  action='store_true',
                        help="Verbose output")

    parser.add_argument("-permute",   "--permute",  action='store_true',
                        help="Permute sequence")

    parser.add_argument("-npermute",  "--npermute",  type=int,
                        help="Number of permutation iterations", default=100)

    parser.add_argument("-writeall",   "--writeall",  action='store_true',
                        help="Write all results")

    parser.add_argument("-seed",   "--seed",  type=int, default=0,
                        help="Random seed")

    args = parser.parse_args()

    # Check if the required arguments are passed to the code
    if args.input is None:
        parser.print_help()
        sys.exit('Input file does not exist!')


if __name__ == "__main__":
    main()
