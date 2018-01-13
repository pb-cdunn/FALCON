import argparse
import collections
import logging
import os
import sys
from .. import io
from .. import bash

LOG = logging.getLogger()

def run(db_prefix, pread_aln, skip_checks, run_jobs_fn, nblock, scattered_fn):
    db_build_done_fn = None
    daligner_scripts = bash.scripts_daligner(run_jobs_fn, db_prefix, db_build_done_fn, nblock=nblock, pread_aln=pread_aln, skip_check=skip_checks)
    basedir = os.path.dirname(os.path.abspath(scattered_fn))
    rootdir = os.path.dirname(os.path.dirname(basedir)) # for now
    jobs = list()
    for job_uid, script in daligner_scripts:
        job_id = 'j_{}'.format(job_uid)
        daligner_settings = dict(db_prefix=db_prefix)

        # Write the scattered inputs.
        daligner_script_fn = '{rootdir}/1-preads_ovl/daligner-scripts/{job_id}/daligner-script.sh'.format(**locals())
        daligner_settings_fn = '{rootdir}/1-preads_ovl/daligner-scripts/{job_id}/settings.json'.format(**locals())
        io.mkdirs(os.path.dirname(daligner_script_fn))
        with open(daligner_script_fn, 'w') as stream:
            stream.write(script)
        io.serialize(daligner_settings_fn, daligner_settings)

        # Record in a job-dict.
        job = dict()
        job['input'] = dict(
                daligner_script = daligner_script_fn,
                daligner_settings = daligner_settings_fn, # not used today, but maybe someday
        )
        job['output'] = dict(
                job_done = '{rootdir}/1-preads_ovl/daligner/{job_id}/daligner.done'.format(**locals()),
        )
        job['params'] = dict(
        )
        job['wildcards'] = {'job_id': job_id} # This should match the wildcard used in the pattern elsewhere.
        jobs.append(job)

    io.serialize(scattered_fn, jobs)

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Prepare for parallel daligner jobs.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--run-jobs-fn',
        help='Input. Result of HPC.daligner.',
    )
    parser.add_argument(
        '--db-prefix', default='preads',
        help='Either preads or raw_reads.',
    )
    parser.add_argument(
        '--skip-checks', default=0,
        help='Skip LAcheck calls after daligner. (0 => do not skip)',
    )
    parser.add_argument(
        '--nblock', default=0,
        help='Number of blocks',
    )
    parser.add_argument(
        '--pread-aln', action='store_true',
        help='Pread alignment mode. (Run daligner_p instead of daligner.)',
    )
    parser.add_argument(
        '--scattered-fn',
        help='Output. JSON list of jobs, where each is a dict of input/output/params/wildcards.',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
