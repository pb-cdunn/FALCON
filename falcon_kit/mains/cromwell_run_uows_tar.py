from __future__ import absolute_import

import argparse
import collections
import glob
import logging
import os
import sys
import pypeflow.do_task
from .. import io

LOG = logging.getLogger()

def dir_from_tar(tar_fn):
    # standard convention for tar-files
    return os.path.splitext(os.path.basename(tar_fn))[0]

def run(tool, uows_tar_fn, nproc):
    cmd = 'tar --strip-components=1 -xvf {}'.format(uows_tar_fn)
    io.syscall(cmd)
    #uows_dn = dir_from_tar(uows_tar_fn)
    uows_dn = '.'
    uows = list(sorted(glob.glob('{}/uow-*'.format(uows_dn))))
    print(uows)
    las_fns = list()
    for uow in uows:
        with io.cd(uow):
            cmd = 'bash -vex uow.sh'
            io.syscall(cmd)
        #las_fns.extend(sorted(glob.glob('{}/*.las'.format(uow))))
    #cmd = 'LAmerge {} {}'.format(
    #    result_fn, ' '.join(las_fns))
    #io.syscall(cmd)
    #io.rm(*las_fns)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Run a bash script once for each unit-of-work, in its own sub-dir. Handle results case-by-case, according to "tool".'
    epilog = '''For now, runs will be in series, since we do not know how many processors we can use.

For tool=daligner, we merge .las files into a single .las
'''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--nproc',
        help='Number of processors to be used.')
    parser.add_argument(
        '--uows-tar-fn',
        help='Input. Tarfile of directories of unit-of-work.')
    parser.add_argument(
        '--tool', default='daligner', choices=['daligner', 'datander'],
        help='The tool for each unit of work. (Currently ignored.)')

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
