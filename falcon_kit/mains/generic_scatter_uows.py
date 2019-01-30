from __future__ import absolute_import

import argparse
import collections
import glob
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def yield_uows(n, all_uows):
    uows_per_chunk = (len(all_uows) + n - 1) / n
    for uow in all_uows:
        yield [uow]


def run(all_uow_list_fn, pattern, nchunks_max):
    all_uows = io.deserialize(all_uow_list_fn)
    n = min(nchunks_max, len(all_uows))
    LOG.info('Num chunks = {} (<= {})'.format(n, nchunks_max))
    all_dn = os.path.abspath(os.path.dirname(all_uow_list_fn))

    for i, uows in enumerate(yield_uows(n, all_uows)):
        key = '{:02d}'.format(i)
        fn = pattern.replace('%', key)
        LOG.info('Writing {} units-of-work to "{}" ({}).'.format(len(uows), fn, key))

        one_dn = os.path.abspath(os.path.dirname(fn))
        rel_dn = os.path.relpath(all_dn, one_dn)
        def fixpath(rel):
            try:
                if not os.path.isabs(rel):
                    return os.path.join('.', os.path.normpath(os.path.join(rel_dn, rel)))
            except Exception:
                # in case of non-string?
                pass
            return rel
        for one_uow in uows:
            if isinstance(one_uow, dict):
                input_dict = one_uow['input']
                for k, v in input_dict.items():
                    input_dict[k] = fixpath(v)

        io.serialize(fn, uows)


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
        '--all-uow-list-fn',
        help='Input. JSON list of all units of work.')
    parser.add_argument(
        '--nchunks-max', type=int,
        help='Input. Maximum number of output files.')
    parser.add_argument(
        '--pattern',
        help='Output. The "%" will be replaced by a zero-padded number. (Probably should be ".json")')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
