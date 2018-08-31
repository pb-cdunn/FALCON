"""
"""
from __future__ import absolute_import
from __future__ import print_function


from future.utils import viewitems
import argparse
import logging
import os
import string
import sys
from ..util import io

LOG = logging.getLogger()


def post_hook(config_fn, db_fn):
    config = io.deserialize(config_fn)
    hook = config.get('LA4Falcon_post')
    if hook:
        LOG.warning('Found LA4Falcon_post in General section of cfg. About to run {!r}...'.format(hook))
        db = os.path.abspath(db_fn)
        parent = os.path.abspath(os.path.dirname(os.getcwd()))
        dbdir = os.path.join(config['LA4Falcon_dbdir'], 'fc-db') + parent
        cmd = string.Template(hook).substitute(DB=db, DBDIR=dbdir)
        io.syscall(cmd)

def run(gathered_fn, db_fn, config_fn, preads_fofn_fn):
    gathered = io.deserialize(gathered_fn)
    d = os.path.abspath(os.path.realpath(os.path.dirname(gathered_fn)))
    def abspath(fn):
        if os.path.isabs(fn):
            return fn # I expect this never to happen though.
        return os.path.join(d, fn)
    fasta_fns = list()
    for desc in gathered:
        fn = abspath(desc['fasta'])
        if 0 == io.filesize(fn):
            LOG.warning('Skipping empty fasta {!r}'.format(fn))
            continue
        fasta_fns.append(fn)
    with open(preads_fofn_fn,  'w') as f:
        for filename in sorted(fasta_fns, key=lambda fn: (os.path.basename(fn), fn)):
            print(filename, file=f)
    post_hook(config_fn, db_fn)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Turn gathered file into FOFN of fasta files.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--gathered-fn',
        help='Input. JSON list of output dicts.')
    parser.add_argument(
        '--db-fn',
        help='Input. Dazzler DB of raw_reads.')
    parser.add_argument(
        '--config-fn',
        help='Input. JSON of relevant configuration (currently from General section of full-prog config).')
    parser.add_argument(
        '--preads-fofn-fn',
        help='Output. FOFN of preads (fasta files).',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
