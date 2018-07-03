from __future__ import absolute_import

from ..pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from .. import run_support
from .. import bash, pype_tasks, snakemake
from ..util.system import (only_these_symlinks, lfs_setstripe_maybe)
from .. import io
from .. import functional
# pylint: disable=no-name-in-module, import-error, fixme, line-too-long
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
                                             makePypeLocalFile, fn, PypeTask)
import argparse
import glob
import json
import logging
import os
import re
import sys
import time


LOG = logging.getLogger(__name__)  # default, for remote tasks

def check_general_config(general_config, input_config_fn):
    if ('pa_daligner_option' not in general_config or
        'ovlp_daligner_option' not in general_config):
        msg = '''Missing options.
We now require both "pa_daligner_option" (stage 0) and "ovlp_daligner_option" (stage 1),
which are automatically passed along to
  HPC.daligner
  HPC.TANmask
  HPC.REPmask

These can provide additional flags:
  pa_HPCdaligner_option
  pa_HPCTANmask_option
  ovlp_HPCdaligner_option
  pa_REPmask_code (-g/-c pairs for 3 iterations, e.g. '1,20;5,15;20,10')
'''
        raise Exception(msg)
    required = ('input_fofn', 'genome_size')
    for name in required:
        assert name in general_config, 'Missing "{}" in {}.'.format(name, input_config_fn)

def main1(prog_name, input_config_fn, logger_config_fn=None):
    global LOG
    LOG = run_support.setup_logger(logger_config_fn)
    lfs_setstripe_maybe(path='.', stripe=12)

    LOG.info('fc_run started with configuration %s', input_config_fn)
    try:
        config = run_support.parse_cfg_file(input_config_fn)
        import json
        dumped = json.dumps(config, indent=2, separators=(',', ': '), sort_keys=True)
        LOG.info('cfg=\n{}'.format(dumped))
    except Exception:
        LOG.exception('Failed to parse config "{}".'.format(input_config_fn))
        raise
    general_config = config['General']
    check_general_config(general_config, input_config_fn)
    input_fofn_plf = makePypeLocalFile(general_config['input_fofn'])
    genome_size = int(general_config['genome_size'])
    squash = True if 0 < genome_size < 1000000 else False
    wf = PypeProcWatcherWorkflow(job_defaults=config['job.defaults'],
                                 squash=squash,
    )
    general_config['ver'] = '100'
    # Store config as JSON, available to many tasks.
    config_fn = './config.json' # must not be in a task-dir
    io.serialize(config_fn, config)
    #with open('foo.snake', 'w') as snakemake_writer:
    with open('/dev/null', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        run(wf, config, rule_writer,
            os.path.abspath(config_fn),
            input_fofn_plf=input_fofn_plf,
            )


def add_bam2dexta_tasks(wf, config, rule_writer,
                        config_fn,
                        input_fofn_plf):
        # run bam2dexta
        bam2dexta_uows_fn = os.path.join(
            rawread_dir, 'bam2dexta-split', 'bam2dexta-uows.json')
        bam2dexta_bash_template_fn = os.path.join(
            rawread_dir, 'bam2dexta-split', 'bash_template.sh')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_BAM2DEXTA_SPLIT_SCRIPT,
            inputs={
                'bam': '/pbi/dept/secondary/testdata/git_sym_cache/synth5k.2016-11-02/synth5k.xml',
            },
            outputs={
                'split': bam2dexta_uows_fn,
                'bash_template': bam2dexta_bash_template_fn,
            },
            parameters={
                'wildcards': 'bam2dexta_id',
            },
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        gathered_fn = os.path.join(rawread_dir, 'bam2dexta-gathered', 'gathered-dexta-files.json')
        gen_parallel_tasks(
            wf, rule_writer,
            bam2dexta_uows_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=bam2dexta_bash_template_fn,
                script='fubar-TODO', #pype_tasks.TASK_DB_TAN_APPLY_SCRIPT, # for snakemake stuff
                inputs={
                    'units_of_work': '0-rawreads/bam2dexta-chunks/{bam2dexta_id}/some-units-of-work.json',
                },
                outputs={
                    'results': '0-rawreads/bam2dexta-runs/{bam2dexta_id}/some-done-files.json',
                },
                parameters={},

            ),
            dist=Dist(NPROC=1, MB=4000, job_dict=config['job.step.da']),
        )

        input_fofn_fn = os.path.join(rawread_dir, 'bam2dexta-combine', 'input.fofn')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_BAM2DEXTA_COMBINE_SCRIPT,
            inputs={
                'gathered': gathered_fn,
            },
            outputs={
                'fofn': input_fofn_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        return gathered_fn, input_fofn_fn

def run(wf, config, rule_writer,
        config_fn,
        input_fofn_plf,
        ):
    """
    Preconditions (for now):
    * LOG
    * run_run_support.logger
    """
    parsed_config = io.deserialize(config_fn)
    if parsed_config != config:
        msg = 'Config from {!r} != passed config'.format(config_fn)
        raise Exception(msg)
    general_config = config['General']
    general_config_fn = os.path.join(os.path.dirname(config_fn), 'General_config.json')
    io.serialize(general_config_fn, general_config) # Some tasks use this.
    rawread_dir = '0-rawreads'
    pread_dir = '1-preads_ovl'
    falcon_asm_dir = '2-asm-falcon'

    for d in (rawread_dir, pread_dir, falcon_asm_dir):
        run_support.make_dirs(d)

    # only matter for parallel jobs
    job_defaults = config['job.defaults']
    #exitOnFailure = bool(job_defaults.get('stop_all_jobs_on_failure', False))
    global default_njobs
    default_njobs = int(job_defaults.get('njobs', 7))
    wf.max_jobs = default_njobs

    assert general_config['input_type'] in (
        'raw', 'preads'), 'Invalid input_type=={!r}'.format(general_config['input_type'])

    parameters = {}

    if general_config['input_type'] == 'raw':
        # Most common workflow: Start with rawreads.

        # import sequences into daligner DB
        # calculate length_cutoff (if specified as -1)
        # split DB
        # run DBdust
        r_db_dust_fn = os.path.join(rawread_dir, 'build', 'raw_reads.db')
        length_cutoff_fn = os.path.join(rawread_dir, 'build', 'length_cutoff')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_BUILD_SCRIPT,
            inputs={
                'config': general_config_fn,
                'input_fofn': fn(input_fofn_plf),
            },
            outputs={
                'length_cutoff': length_cutoff_fn,
                'db': r_db_dust_fn,
                # Also .raw_reads.*, of course. And dust track.
            },
            parameters=dict(
            ),
            rule_writer=rule_writer,
            dist=Dist(NPROC=1),
        ))

        # run TANmask
        tan_uows_fn = os.path.join(
            rawread_dir, 'tan-split', 'tan-uows.json')
        tan_bash_template_fn = os.path.join(
            rawread_dir, 'tan-split', 'bash_template.sh')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_TAN_SPLIT_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': r_db_dust_fn,
            },
            outputs={
                'split': tan_uows_fn,
                'bash_template': tan_bash_template_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(NPROC=1),
        ))

        gathered_fn = os.path.join(rawread_dir, 'tan-gathered', 'gathered-done-files.json')
        gen_parallel_tasks(
            wf, rule_writer,
            tan_uows_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=tan_bash_template_fn,
                script='fubar-TODO', #pype_tasks.TASK_DB_TAN_APPLY_SCRIPT, # for snakemake stuff
                inputs={
                    'units_of_work': '0-rawreads/tan-chunks/{tan0_id}/some-units-of-work.json',
                },
                outputs={
                    #'job_done': '0-rawreads/{dal0_id}/daligner.done',
                    'results': '0-rawreads/tan-runs/{tan0_id}/some-done-files.json',
                },
                parameters={},

            ),
            dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.da']),
        )

        r_db_tan_fn = os.path.join(rawread_dir, 'tan-combine', 'raw_reads.db')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_TAN_COMBINE_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': r_db_dust_fn,
                'gathered': gathered_fn,
            },
            outputs={
                'new_db': r_db_tan_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        #### HPC.REPmask/daligner/LAmerge
        codes = functional.parse_REPmask_code(general_config['pa_REPmask_code'])
        LOG.info('Parsed pa_REPmask_code (repa,repb,repc): {!r}'.format(codes))

        ### REPmask tasks (a, b, c)
        letter = 'a'
        group_size, coverage_limit = codes[0]
        i_db_fn = r_db_tan_fn
        o_db_fn = add_rep_tasks(wf, rule_writer, rawread_dir, config, general_config,
                general_config_fn, i_db_fn, length_cutoff_fn,
                letter, group_size, coverage_limit)
        letter = 'b'
        group_size, coverage_limit = codes[1]
        i_db_fn = o_db_fn
        o_db_fn = add_rep_tasks(wf, rule_writer, rawread_dir, config, general_config,
                general_config_fn, i_db_fn, length_cutoff_fn,
                letter, group_size, coverage_limit)
        letter = 'c'
        group_size, coverage_limit = codes[2]
        i_db_fn = o_db_fn
        o_db_fn = add_rep_tasks(wf, rule_writer, rawread_dir, config, general_config,
                general_config_fn, i_db_fn, length_cutoff_fn,
                letter, group_size, coverage_limit)
        r_db_rep_fn = o_db_fn

        #### basic daligner/LAmerge
        p_id2las_fn = os.path.join(rawread_dir, 'las-merge-combine', 'p_id2las.json')
        las_fofn_fn = os.path.join(rawread_dir, 'las-merge-combine', 'las_fofn.json')

        add_daligner_and_merge_tasks(
            wf, rule_writer,
            general_config, config['job.step.da'], config['job.step.la'],
            rawread_dir,
            general_config_fn, r_db_rep_fn,
            length_cutoff_fn,
            p_id2las_fn, las_fofn_fn,
            daligner_wildcard='dal0_id',
            lamerge_wildcard='mer0_id',
            daligner_params={},
            db_prefix='raw_reads', # TODO: Infer
            daligner_split_script=pype_tasks.TASK_DB_DALIGNER_SPLIT_SCRIPT,
        )
        ####

        if general_config['target'] == 'overlapping':
            sys.exit(0)

        # Produce new FOFN of preads fasta, based on consensus of overlaps.
        wf.max_jobs = config['job.step.cns'].get('njobs', default_njobs)

        split_fn = os.path.join(
            rawread_dir, 'cns-split', 'split.json')
        bash_template_fn = os.path.join(
            rawread_dir, 'cns-split', 'consensus-bash-template.sh')
        params = dict(parameters)
        params['wildcards'] = 'cns0_id,cns0_id2'
        wf.addTask(gen_task(
            script=pype_tasks.TASK_CONSENSUS_SPLIT_SCRIPT,
            inputs={
                'p_id2las': p_id2las_fn,
                'raw_reads_db': r_db_rep_fn,
                'length_cutoff': length_cutoff_fn,
                'config': general_config_fn,
            },
            outputs={
                'split': split_fn,
                'bash_template': bash_template_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        gathered_fn = os.path.join(rawread_dir, 'cns-gather', 'gathered.json')
        gen_parallel_tasks(
            wf, rule_writer,
            split_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=bash_template_fn,
                script=pype_tasks.TASK_CONSENSUS_TASK_SCRIPT, # for snakemake only
                inputs = {
                    #'las': '0-rawreads/cns-split/{cns0_id}/merged.{cns0_id2}.las',
                    #'db': r_db_rep_fn,
                    #'length_cutoff': length_cutoff_fn,
                    #'config': general_config_fn,
                    'units_of_work': '0-rawreads/cns-chunks/{cns0_id}/some-units-of-work.json',
                },
                outputs = {
                    #'fasta': '0-rawreads/consensus/{cns0_id}/consensus.{cns0_id2}.fasta',
                    'results': '0-rawreads/cns-runs/{cns0_id}/some-done-files.json',
                },
                parameters={},
            ),
            dist=Dist(NPROC=6, job_dict=config['job.step.cns']),
        )
        preads_fofn_fn = os.path.join(rawread_dir, 'preads', 'input_preads.fofn')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_CONSENSUS_GATHER_SCRIPT,
            inputs={
                'gathered': gathered_fn,
            },
            outputs={
                'preads_fofn': preads_fofn_fn,
            },
            parameters=parameters, #{},
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        rdir = os.path.join(rawread_dir, 'report')
        pre_assembly_report_fn = os.path.join(rdir, 'pre_assembly_stats.json')
        params = dict(parameters)
        params['length_cutoff_user'] = general_config['length_cutoff']
        params['genome_length'] = general_config['genome_size'] # note different name; historical
        wf.addTask(gen_task(
            script=pype_tasks.TASK_REPORT_PRE_ASSEMBLY_SCRIPT,
            inputs={'length_cutoff': length_cutoff_fn,
                    'raw_reads_db': r_db_rep_fn,
                    'preads_fofn': preads_fofn_fn,
                    'config': general_config_fn,
            },
            outputs={'pre_assembly_report': pre_assembly_report_fn,
            },
            parameters=params,
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        if general_config['target'] == 'pre-assembly':
            wf.refreshTargets()
            LOG.info('Quitting after stage-0 for General.target=pre-assembly')
            return

    # build pread database
    if general_config['input_type'] == 'preads':
        LOG.info('General.input_type=preads, so we skip stage 0-rawreads.')
        preads_fofn_fn = general_config['input_fofn']
        assert os.path.exists(preads_fofn_fn), '{!r} does not exist.'.format(preads_fofn_fn)

    pdb_build_done = os.path.join(pread_dir, 'pdb_build_done')
    run_jobs_fn = os.path.join(pread_dir, 'run_jobs.sh')
    preads_db_fn = os.path.join(pread_dir, 'build', 'preads.db')
    length_cutoff_pr_fn = os.path.join(pread_dir, 'build', 'length_cutoff')

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_BUILD_SCRIPT,
        inputs={
            'config': general_config_fn,
            'input_fofn': preads_fofn_fn,
        },
        outputs={
            'length_cutoff': length_cutoff_pr_fn,
            'db': preads_db_fn,
            # Also .preads.*, of course.
        },
        parameters=dict(
        ),
        rule_writer=rule_writer,
        dist=Dist(NPROC=1),
    ))

    ####
    p_id2las_fn = os.path.join(pread_dir, 'las-merge-combine', 'block2las.json')
    las_fofn_fn = os.path.join(pread_dir, 'las-merge-combine', 'las_fofn.json')

    add_daligner_and_merge_tasks(
        wf, rule_writer,
        general_config, config['job.step.pda'], config['job.step.pla'],
        pread_dir,
        general_config_fn, preads_db_fn, # no tan-mask for preads
        length_cutoff_pr_fn,
        p_id2las_fn, las_fofn_fn,
        daligner_wildcard='dal1_id',
        lamerge_wildcard='mer1_id',
        daligner_params={},
        db_prefix='preads', # TODO: Infer
        daligner_split_script=pype_tasks.TASK_DB_DALIGNER_SPLIT_SCRIPT,
    )
    ####

    wf.max_jobs = config['job.step.asm'].get('njobs', default_njobs)
    db2falcon_dir = os.path.join(pread_dir, 'db2falcon')
    db2falcon_done_fn = os.path.join(db2falcon_dir, 'db2falcon_done')
    preads4falcon_fn = os.path.join(db2falcon_dir, 'preads4falcon.fasta')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_RUN_DB_TO_FALCON_SCRIPT,
        inputs={'p_id2las': p_id2las_fn,
                'preads_db': preads_db_fn,
                },
        outputs={'job_done': db2falcon_done_fn,
                 'preads4falcon': preads4falcon_fn,
                 },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(NPROC=4, job_dict=config['job.step.asm']),
    ))

    falcon_asm_done_fn = os.path.join(falcon_asm_dir, 'falcon_asm_done')
    for key in ('overlap_filtering_setting', 'length_cutoff_pr', 'fc_ovlp_to_graph_option'):
        parameters[key] = general_config[key]
    wf.addTask(gen_task(
        script=pype_tasks.TASK_RUN_FALCON_ASM_SCRIPT,
        inputs={'db2falcon_done': db2falcon_done_fn, 'db_file': preads_db_fn,
                'preads4falcon_fasta': preads4falcon_fn,
                'las_fofn': las_fofn_fn,
                'config': general_config_fn,
                },
        outputs={'falcon_asm_done': falcon_asm_done_fn},
        parameters=parameters,
        rule_writer=rule_writer,
        dist=Dist(NPROC=4, job_dict=config['job.step.asm']),
    ))
    wf.refreshTargets()

    with io.cd('0-rawreads'):
        # for backwards-compatibility
        io.symlink('las-merge-combine', 'las-gather')

    #return falcon_asm_done
def add_daligner_and_merge_tasks(
        wf, rule_writer,
        general_config, daligner_job_config, merge_job_config,
        super_dir,
        general_config_fn, db_fn,
        length_cutoff_fn, # not always needed (refactor later)
        p_id2las_fn, las_fofn_fn,
        daligner_wildcard, #='dal0_id',
        lamerge_wildcard, #='mer0_id',
        daligner_params=dict(),
        db_prefix='raw_reads',
        daligner_split_script=pype_tasks.TASK_DB_DALIGNER_SPLIT_SCRIPT,
    ):
    """
    Results:
      p_id2las_fn, las_fofn_fn
    """
    max_jobs = wf.max_jobs
    parameters = dict()

    # run daligner
    wf.max_jobs = daligner_job_config.get('njobs', default_njobs)
    daligner_all_units_fn = os.path.join(
        super_dir, 'daligner-split', 'all-units-of-work.json')
    daligner_bash_template_fn = os.path.join(
        super_dir, 'daligner-split', 'daligner_bash_template.sh')
    params = dict(daligner_params)
    params['skip_checks'] = int(general_config.get('skip_checks', 0))
    params['wildcards'] = daligner_wildcard
    wf.addTask(gen_task(
        script=daligner_split_script,
        inputs={
            'config': general_config_fn,
            'db': db_fn,
            'length_cutoff': length_cutoff_fn,
        },
        outputs={
            'split': daligner_all_units_fn,
            'bash_template': daligner_bash_template_fn
        },
        parameters=params,
        rule_writer=rule_writer,
        dist=Dist(local=True, NPROC=4), # really, NPROC=1, but we need to know the max
    ))

    gathered_fn = os.path.join(super_dir, 'daligner-gathered', 'gathered-done-files.json')
    gen_parallel_tasks(
        wf, rule_writer,
        daligner_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=daligner_bash_template_fn,
            script=pype_tasks.TASK_DB_DALIGNER_APPLY_SCRIPT, # for snakemake stuff
            inputs={
                'units_of_work': os.path.join(super_dir, 'daligner-chunks/{%s}/some-units-of-work.json'%daligner_wildcard),
            },
            outputs={
                'results': os.path.join(super_dir, 'daligner-runs/{%s}/some-done-files.json'%daligner_wildcard),
            },
            parameters={},
        ),
        dist=Dist(NPROC=4, MB=4000, job_dict=daligner_job_config),
    )

    gathered_las_fn = os.path.join(super_dir, 'daligner-combine', 'gathered-las.json')
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_DALIGNER_COMBINE_SCRIPT,
        inputs={
            'config': general_config_fn,
            'db': db_fn,
            'gathered': gathered_fn,
        },
        outputs={
            'las_paths': gathered_las_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        #dist=Dist(NPROC=1, MB=4000, job_dict=daligner_job_config)
        dist=Dist(local=True),
    ))

    # Merge .las files.
    wf.max_jobs = merge_job_config.get('njobs', default_njobs)
    las_merge_all_units_fn = os.path.join(super_dir, 'las-merge-split', 'all-units-of-work.json')
    bash_template_fn = os.path.join(super_dir, 'las-merge-split', 'las-merge-bash-template.sh')
    params = dict(parameters)
    params['db_prefix'] = db_prefix
    params['wildcards'] = lamerge_wildcard
    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_LAMERGE_SPLIT_SCRIPT,
        inputs={
            'config': general_config_fn,
            'las_paths': gathered_las_fn,
        },
        outputs={
            'split': las_merge_all_units_fn,
            'bash_template': bash_template_fn,
        },
        parameters=params,
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    gathered_fn = os.path.join(super_dir, 'las-merge-gathered', 'gathered.json')
    gen_parallel_tasks(
        wf, rule_writer,
        las_merge_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=bash_template_fn,
            script=pype_tasks.TASK_DB_LAMERGE_APPLY_SCRIPT, # for snakemake
            inputs={
                'units_of_work': os.path.join(super_dir, 'las-merge-chunks/{%s}/some-units-of-work.json'%lamerge_wildcard),
            },
            outputs={
                'results': os.path.join(super_dir, 'las-merge-runs/{%s}/some-las-paths.json'%lamerge_wildcard),
            },
            parameters={},
        ),
        dist=Dist(NPROC=1, job_dict=merge_job_config),
    )

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DB_LAMERGE_COMBINE_SCRIPT,
        inputs={
            'config': general_config_fn,
            'gathered': gathered_fn,
        },
        outputs={
            'block2las': p_id2las_fn,
            'las_paths': las_fofn_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    wf.max_jobs = max_jobs


def add_rep_tasks(
        wf, rule_writer,
        rawread_dir, config, general_config,
        general_config_fn, i_db_fn, length_cutoff_fn,
        letter, group_size, coverage_limit,
        ):
        """
        Add daligner/lamerge/REPmask parallel tasks for one iteration of repeat-masking.
        TODO: Make the tasks no-ops if the codes are zero (or something like that).
        """
        name = 'rep{}'.format(letter)
        rep_dir = os.path.join(rawread_dir, name)
        o_db_rep_fn = os.path.join(rep_dir, 'rep-combine', 'raw_reads.db')

        p_id2las_fn = os.path.join(rep_dir, 'las-merge-combine', 'p_id2las.json')
        las_fofn_fn = os.path.join(rep_dir, 'las-merge-combine', 'las_fofn.json')

        rep_daligner_params = dict(
            group_size=group_size, coverage_limit=coverage_limit,
        )
        add_daligner_and_merge_tasks(
            wf, rule_writer,
            general_config, config['job.step.da'], config['job.step.la'],
            rep_dir,
            general_config_fn, i_db_fn,
            length_cutoff_fn,
            p_id2las_fn, las_fofn_fn,
            daligner_wildcard='dal0{}_id'.format(letter),
            lamerge_wildcard='mer0{}_id'.format(letter),
            daligner_params=rep_daligner_params,
            db_prefix='raw_reads', # TODO: Infer
            daligner_split_script=pype_tasks.TASK_DB_REP_DALIGNER_SPLIT_SCRIPT,
        )

        ### REPmask
        # rep-split
        # We assume that daligner/LAmerge have already run.
        # Instead of using the REP.mask calls from rep-jobs.05.MASK,
        # we construct our own.
        rep_uows_fn = os.path.join(
            rep_dir, 'rep-split', 'rep-uows.json')
        rep_bash_template_fn = os.path.join(
            rep_dir, 'rep-split', 'bash_template.sh')
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_REP_SPLIT_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': i_db_fn,
                'las_paths': las_fofn_fn,
            },
            outputs={
                'split': rep_uows_fn,
                'bash_template': rep_bash_template_fn,
            },
            parameters={
                'group_size': group_size,
                'coverage_limit': coverage_limit,
                'wildcards': '{}_id'.format(name),
            },
            rule_writer=rule_writer,
            dist=Dist(NPROC=1),
        ))

        # rep-apply
        gathered_fn = os.path.join(rep_dir, 'rep-gathered', 'gathered-done-files.json')
        gen_parallel_tasks(
            wf, rule_writer,
            rep_uows_fn, gathered_fn,
            run_dict=dict(
                bash_template_fn=rep_bash_template_fn,
                script='fubar-TODO', #pype_tasks.TASK_DB_REP_APPLY_SCRIPT, # for snakemake stuff
                inputs={
                    'units_of_work': '0-rawreads/%(name)s/rep-chunks/{%(name)s_id}/some-units-of-work.json'%locals(),
                },
                outputs={
                    'results': '0-rawreads/%(name)s/rep-runs/{%(name)s_id}/some-done-files.json'%locals(),
                },
                parameters={},

            ),
            dist=Dist(NPROC=4, MB=4000, job_dict=config['job.step.da']),
        )

        # rep-combine
        wf.addTask(gen_task(
            script=pype_tasks.TASK_DB_REP_COMBINE_SCRIPT,
            inputs={
                'config': general_config_fn,
                'db': i_db_fn,
                'gathered': gathered_fn,
            },
            outputs={
                'new_db': o_db_rep_fn,
            },
            parameters={
                'group_size': group_size,
            },
            rule_writer=rule_writer,
            dist=Dist(local=True),
        ))

        return o_db_rep_fn


def main(argv=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('config',
                        help='.cfg/.ini/.json')
    parser.add_argument('logger',
                        nargs='?',
                        help='(Optional)JSON config for standard Python logging module')
    args = parser.parse_args(argv[1:])
    main1(argv[0], args.config, args.logger)


if __name__ == '__main__':
    main()
