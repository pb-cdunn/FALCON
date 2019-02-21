from __future__ import absolute_import

import argparse
import glob
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def rename(fn, prefix):
    dn, bn = os.path.split(fn)
    nfn = os.path.join(dn, prefix + bn)
    cmd = 'mv -f {} {}'.format(fn, nfn)
    LOG.info('cmd ={!r}'.format(cmd))
    io.syscall(cmd)

def run(pattern, prefix):
    for fn in glob.iglob(pattern):
        rename(fn, prefix)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Find all files matching "pattern" under CWD, and prefix with "prefix".'
    epilog = 'We do this because Cromwell does not properly glob dot-files.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--pattern',
        help='Find all files matching this, including in subdirs.')
    parser.add_argument(
        '--prefix', default='dot',
        help='Rename the matching files to have this prefix. E.g. ".foo" becomes "PREFIX.foo".')

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
