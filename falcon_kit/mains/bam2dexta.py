from __future__ import absolute_import

import argparse
import collections
import glob
import logging
import os
import re
import sys
import time
from ..util.io import yield_validated_fns
from .. import io, functional
from .. import(
        bash,  # for write_sub_script
        pype_tasks,  # for TASKS
)

LOG = logging.getLogger()
WAIT = 20 # seconds


def bam2dexta_split(config_fn, bam_subreadset_fn, wildcards, split_fn, bash_template_fn):
    assert bam_subreadset_fn.endswith('.xml')
    with open(bash_template_fn, 'w') as stream:
        stream.write(pype_tasks.TASK_BAM2DEXTA_APPLY_SCRIPT)

    split_dataset_prefix = os.path.join(os.getcwd(), 'split') # TODO: Test this as relative sub-dir.

    from ..util import dataset_split # introduces pbcore dependency
    bam_paths = dataset_split.split_dataset(bam_subreadset_fn, split_dataset_prefix)

    jobs = list()
    for i, bam_fn in enumerate(bam_paths):
        job_id = 'b_{:03d}'.format(i)

        # Write the las files for this job.
        #input_dir = os.path.join('bam2dexta-scripts', job_id)
        #bam_paths_fn = os.path.join('.', input_dir, 'bam-paths.json')
        #io.mkdirs(input_dir)
        #io.serialize(bam_paths_fn, bam_paths)

        # Record in a job-dict.
        dexta_fn = 'subreads.{}.dexta'.format(job_id)
        job = dict()
        job['input'] = dict(
                config=config_fn,
                bam=bam_fn,
        )
        job['output'] = dict(
                dexta=dexta_fn
        )
        job['params'] = dict(
        )
        job['wildcards'] = {wildcards: job_id}
        jobs.append(job)
    io.serialize(split_fn, jobs)

def bam2dexta_apply(bam_fn, dexta_fn):
    """Given a bam subread DataSet, write a .dexta file.
    """
    io.rm_force(dexta_fn)
    # This command is not quite right yet. TODO
    cmd = 'rm -f {dexta_fn}; bam2fasta -u -o foo {bam_fn}; dexta foo.fasta; time mv -f foo.dexta {dexta_fn}'.format(
            **locals())
    # Note: If 'dexta' fails, the script will error. So we might still have an empty foo.dexta, but
    # we will not have moved it to {dexta_fn}.
    io.syscall(cmd)

def bam2dexta_combine(gathered_fn, dexta_fofn_fn):
    gathered = io.deserialize(gathered_fn)
    d = os.path.abspath(os.path.realpath(os.path.dirname(gathered_fn)))
    def abspath(fn):
        if os.path.isabs(fn):
            return fn # I expect this never to happen though.
        return os.path.join(d, fn)
    dexta_fns = list()
    for job_output in gathered:
        assert len(job_output) == 1, 'len(job_output) == {} != 1'.format(len(job_output))
        for fn in job_output.values():
            abs_fn = abspath(fn)
            dexta_fns.append(abs_fn)
    dexta_paths = list()
    for dexta_fn in sorted(dexta_fns):
        if not os.path.exists(dexta_fn):
            msg = 'Did not find {!r}. Waiting {} seconds.'.format(dexta_fn, WAIT)
            LOG.info(msg)
            time.sleep(WAIT)
            if not os.path.exists(dexta_fn):
                msg = 'Did not find {!r}, even after waiting {} seconds. Maybe retry later?'.format(dexta_fn, WAIT)
                raise Exception(msg)
        dexta_paths.append(dexta_fn)

    # Serialize result.
    #io.serialize(dexta_paths_fn, sorted(dexta_paths))
    with open(dexta_fofn_fn, 'w') as stream:
        stream.write('\n'.join(dexta_paths))
        stream.write('\n')


def setup_logging(log_level):
    hdlr = logging.StreamHandler(sys.stderr)
    hdlr.setLevel(log_level)
    hdlr.setFormatter(logging.Formatter('[%(levelname)s]%(message)s'))
    LOG.addHandler(hdlr)
    LOG.setLevel(logging.NOTSET)
    LOG.info('Log-level: {}'.format(log_level))

def cmd_split(args):
    bam2dexta_split(
            args.config_fn, args.bam_subreadset_fn,
            args.wildcards,
            args.split_fn, args.bash_template_fn,
    )
def cmd_apply(args):
    bam2dexta_apply(args.bam_fn, args.dexta_fn)
def cmd_combine(args):
    bam2dexta_combine(args.gathered_fn, args.dexta_fofn_fn)

def get_ours(config_fn, db_fn):
    ours = dict()
    config = io.deserialize(config_fn)
    LOG.info('config({!r}):\n{}'.format(config_fn, config))
    LOG.info('our subset of config:\n{}'.format(ours))
    return ours

def add_split_arguments(parser):
    parser.add_argument(
        '--wildcards', default='bam2dexta0_id',
        help='Comma-separated string of keys to be subtituted into output paths for each job, if any. (Helps with snakemake and pypeflow; not needed in pbsmrtpipe, since outputs are pre-determined.)',
    )
    parser.add_argument(
        '--bam-subreadset-fn',
        help='input. Dataset (.xml) of bam files of subreads.'
    )
    parser.add_argument(
        '--split-fn', default='bam2dexta-uows.json',
        help='output. Units-of-work for bam2fasta/dexta.',
    )
    parser.add_argument(
        '--bash-template-fn', default='bash-template.sh',
        help='output. Script to apply later.',
    )
def add_apply_arguments(parser):
    parser.add_argument(
        '--bam-fn', required=True,
        help='input. bam or dataset')
    parser.add_argument(
        '--dexta-fn', required=True,
        help='output. The dazzler (Gene Myers) dexta-file.',
    )
def add_combine_arguments(parser):
    parser.add_argument(
        '--gathered-fn', required=True,
        help='input. List of sentinels. Produced by gen_parallel_tasks() gathering. The .las files are next to these.',
    )
    parser.add_argument(
        '--dexta-fofn-fn', required=True,
        help='output. FOFN of dexta paths.')

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_args(argv):
    description = 'Efficiently generate .dexta from BAM or subread datasets.'
    epilog = 'For more details on .dexta, see https://dazzlerblog.wordpress.com/command-guides/dextractor-command-guide/'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--log-level', default='INFO',
        help='Python logging level.',
    )
    parser.add_argument(
        '--nproc', type=int, default=0,
        help='ignored for now, but non-zero will mean "No more than this."',
    )
    parser.add_argument(
        '--config-fn', required=True,
        help='Input. JSON of user-configuration. (This is probably the [General] section.)',
    )

    help_split = 'get each bam-file (or subread dataset file)'
    help_apply = 'run bam2fasta and dexta as a unit-of-work'
    help_combine = 'generate a file of .dexta files'

    subparsers = parser.add_subparsers(help='sub-command help')

    parser_split = subparsers.add_parser('split',
            formatter_class=HelpF,
            description=help_split,
            epilog='',
            help=help_split)
    add_split_arguments(parser_split)
    parser_split.set_defaults(func=cmd_split)

    parser_apply = subparsers.add_parser('apply',
            formatter_class=HelpF,
            description=help_apply,
            epilog='',
            help=help_apply)
    add_apply_arguments(parser_apply)
    parser_apply.set_defaults(func=cmd_apply)

    parser_combine = subparsers.add_parser('combine',
            formatter_class=HelpF,
            description=help_combine,
            epilog='',
            help=help_combine)
    add_combine_arguments(parser_combine)
    parser_combine.set_defaults(func=cmd_combine)

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    setup_logging(args.log_level)
    args.func(args)


if __name__ == '__main__':  # pragma: no cover
    main()
