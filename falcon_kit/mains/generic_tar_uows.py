from __future__ import absolute_import, division

import argparse
import collections
import glob
import logging
import os
import sys
import pypeflow.do_task
from .. import io

LOG = logging.getLogger()

def tar_uows(fn, uows):
    # Operate in a subdir. (Named, so not thread-safe.)
    subdir = os.path.splitext(fn)[0]
    io.mkdirs(subdir) # permissions?
    with io.cd(subdir):
        # We could include other files here, or at least symlinks, but not today.
        # Soon, we will construct the uow-subdirs here, but we must consider clobbering.
        io.serialize('some-units-of-work.json', uows)
    cmd = 'tar -cf {} {}'.format(fn, subdir)
    io.syscall(cmd)
    io.rmdirs(subdir)

def yield_uows(n, all_uows):
    uows_per_chunk = (len(all_uows) + n - 1) / n
    for uow in all_uows:
        yield [uow]

def run(all_uow_list_fn, pattern, nchunks_max):
    all_uows = io.deserialize(all_uow_list_fn)
    n = min(nchunks_max, len(all_uows))
    LOG.info('Num chunks = {} (<= {})'.format(n, nchunks_max))
    for i, uows in enumerate(yield_uows(n, all_uows)):
        key = '{:02d}'.format(i)
        fn = pattern.replace('%', key)
        LOG.info('Writing {} units-of-work to "{}" ({}).'.format(len(uows), fn, key))
        tar_uows(fn, uows)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Split a JSON list of units-of-work into up to N files ("chunks"), still as lists of units-of-work.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--all-uow-list-fn',
        help='Input. JSON list of units of work.')
    parser.add_argument(
        '--nchunks_max', type=int,
        help='Input. Maximum number of output files.')
    parser.add_argument(
        '--pattern',
        help='Output. The "%" will be replace by a zero-padded number. (These will be a tar-files, so it should probably end in ".tar".')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
