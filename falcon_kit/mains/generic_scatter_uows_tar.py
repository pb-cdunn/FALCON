from __future__ import absolute_import, division

import argparse
import collections
import glob
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def yield_uows(n, all_uows):
    """
    >>> list(yield_uows(2, [0,1,2,3,4]))
    [[0, 1, 2], [3, 4]]
    """
    # yield exactly n sublists.
    total = len(all_uows)
    remaining = total
    it = iter(all_uows)
    while n and remaining:
        taken = (remaining + n - 1) // n
        to_yield = [next(it) for _ in range(taken)]
        yield to_yield
        remaining -= taken
        n -= 1


def move_into_tar(dn, fns):
    # Create directory 'dn'.
    # Move files (or dir-trees) into directory 'dn', and tar it.
    # By convention, for tar-file "foo.tar", we first move everything into a directory named "foo".
    io.mkdirs(dn)
    for fn in fns:
        cmd = 'mv {} {}'.format(fn, dn)
        io.syscall(cmd)
    tar_fn = '{}.tar'.format(dn)
    #with tarfile.TarFile(tar_fn, 'w', dereference=False, ignore_zeros=True, errorlevel=2) as tf:
    #    tf.add(dn)
    cmd = 'tar cvf {} {}'.format(tar_fn, dn)
    io.syscall(cmd)
    io.rmdirs(dn)


def dir_from_tar(tar_fn):
    return os.path.splitext(os.path.basename(tar_fn))[0]


def run(all_uows_tar_fn, pattern, nchunks_max):
    cmd = 'tar -xvf {}'.format(all_uows_tar_fn)
    io.syscall(cmd)
    all_uows_dn = dir_from_tar(all_uows_tar_fn)
    all_uows = list(sorted(glob.glob('{}/uow-*'.format(all_uows_dn))))
    n = min(nchunks_max, len(all_uows))
    LOG.info('Num chunks = {} (<= {})'.format(n, nchunks_max))

    for i, uows in enumerate(yield_uows(n, all_uows)):
        key = '{:02d}'.format(i)
        fn = pattern.replace('%', key)
        LOG.info('Writing {} units-of-work to "{}" ({}).'.format(len(uows), fn, key))
        dn = dir_from_tar(fn)
        move_into_tar(dn, uows)
    cmd = 'rmdir {}'.format(all_uows_dn)
    io.syscall(cmd)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Scatter a single unit-of-work from many units-of-work.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--all-uows-tar-fn',
        help='Input. Tarfile of all units of work directories.')
    parser.add_argument(
        '--nchunks-max', type=int,
        help='Input. Maximum number of output files.')
    parser.add_argument(
        '--pattern',
        help='Output. The "%" will be replaced by a zero-padded number. (Probably should be ".tar")')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
