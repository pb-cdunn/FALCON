"""Given a FOFN, write JSON list.

Cromwell write_json() does not work as we would expect.

https://github.com/broadinstitute/cromwell/issues/4625

So we use write_lines() instead.

Then, this little program can convert those lines into JSON.
"""
from __future__ import absolute_import

import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(lines_fn, json_fn):
    with open(lines_fn) as sin:
        fns = [line.strip() for line in sin]
    io.serialize(json_fn, fns)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Symlink into current directory. This helps keep command-lines short later.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--lines-fn',
        help='Input. Result of WDL write_lines().')
    parser.add_argument(
        '--json-fn',
        help='Output. Should have been result of WDL write_json().')

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
