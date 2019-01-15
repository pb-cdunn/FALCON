import falcon_kit.mains.fasta_subsample as mod
import functools
import helpers
import pytest
import os
from falcon_kit.FastaReader import open_fasta_reader
import cStringIO
import json
import falcon_kit.mains.fasta_filter as fasta_filter

import random

RANDOM_SEED = 12345

tests_strategy_random = [
    {   # Test empty input.
        'input_data': [],
        'coverage': 1,
        'genome_size': 1,
        'strategy': mod.STRATEGY_RANDOM,

        'exp_whitelist': '[]',
        'exp_all': '[]',
        'exp_stats': ''
    },
    {   # Test on an input with 2 subreads.
       'input_data': ["""\
>synthetic/1/500_1000
GATTACAGATTACAGATTACAGATTACA
>synthetic/2/0_4
ACTG
"""],
        'coverage': 1,
        'genome_size': 2,
        'strategy': mod.STRATEGY_RANDOM,

        'exp_whitelist': '["synthetic/2"]',
        'exp_all': '[["synthetic/1", 28], ["synthetic/2", 4]]',
        'exp_stats': ''
    },
    {   # Test a more complex input of several files, and multiple subreads per ZMW.
       'input_data': ["""\
>synthetic/1/0_5
ACGTAACGTA
>synthetic/1/10_14
ACGTACGTA
>synthetic/1/14_15
ACGTA
""",
"""\
>synthetic/2/0_500
GATTACAGATTACA
>synthetic/2/0_500
GATTACAGATTACAGATTACAGATTACA
>synthetic/2/0_500
GATTACAGATTACAGATTACA
>synthetic/2/500_1000
GATTACAGATTACAGATTACAGATTACA
>synthetic/2/0_500
GATTACA
""",
"""\
>synthetic/3/0_8
ACGTACGTA
>synthetic/4/0_8
ACGTACGTAAAAAAA
>synthetic/5/0_8
ACGTACGTAAAAAAAAAAAA
"""],
        'coverage': 3,
        'genome_size': 10,
        'strategy': mod.STRATEGY_RANDOM,

        'exp_whitelist': '["synthetic/2", "synthetic/4"]',
        'exp_all': '[["synthetic/1", 9], ["synthetic/2", 21], ["synthetic/3", 9], ["synthetic/4", 15], ["synthetic/5", 20]]',
        'exp_stats': ''
    },
]

tests_strategy_longest = [
    {
        'input_data': tests_strategy_random[0]['input_data'],
        'coverage': tests_strategy_random[0]['coverage'],
        'genome_size': tests_strategy_random[0]['genome_size'],
        'strategy': mod.STRATEGY_LONGEST,
        'exp_whitelist': '[]',
        'exp_all': '[]',
        'exp_stats': ''
    },
    {
        'input_data': tests_strategy_random[1]['input_data'],
        'coverage': tests_strategy_random[1]['coverage'],
        'genome_size': tests_strategy_random[1]['genome_size'],
        'strategy': mod.STRATEGY_LONGEST,
        'exp_whitelist': '["synthetic/1"]',
        'exp_all': '[["synthetic/1", 28], ["synthetic/2", 4]]',
        'exp_stats': ''
    },
    {
        'input_data': tests_strategy_random[2]['input_data'],
        'coverage': tests_strategy_random[2]['coverage'],
        'genome_size': tests_strategy_random[2]['genome_size'],
        'strategy': mod.STRATEGY_LONGEST,
        'exp_whitelist': '["synthetic/2", "synthetic/5"]',
        'exp_all': '[["synthetic/1", 9], ["synthetic/2", 21], ["synthetic/3", 9], ["synthetic/4", 15], ["synthetic/5", 20]]',
        'exp_stats': ''
    },
]

def helper_create_input_files(tmpdir, input_data):
    # Create the input files from the given sequences.
    input_files = []
    for i, input_seqs in enumerate(input_data):
        fn = tmpdir.join('in-{}'.format(i))
        fn.write(input_seqs)
        input_files.append(str(fn))
    fn = tmpdir.join('input.fofn')
    fn.write('\n'.join(input_files))
    return input_files, str(fn)

def check_run(tmpdir, input_data, coverage, genome_size, strategy, exp_whitelist, exp_all, exp_stats):
    strategy_func = mod.get_strategy_func(strategy)
    random.seed(RANDOM_SEED)

    input_files, input_fofn_fn = helper_create_input_files(tmpdir, input_data)
    zmws_whitelist, zmws_all, stats_dict = mod.run(fasta_filter.yield_zmwtuple(mod.yield_record(input_files), None, False), coverage, genome_size, strategy_func)

    assert sorted(json.loads(exp_whitelist)) == sorted(zmws_whitelist)
    assert sorted([tuple(val) for val in json.loads(exp_all)]) == sorted(zmws_all)

def check_main(tmpdir, input_data, coverage, genome_size, strategy, exp_whitelist, exp_all, exp_stats):
    out_prefix = tmpdir.join('out')
    random.seed(RANDOM_SEED)

    input_files, input_fofn_fn = helper_create_input_files(tmpdir, input_data)

    argv = ['prog', '--strategy', strategy,
            '--coverage', str(coverage),
            '--genome-size', str(genome_size),
            '--random-seed', str(RANDOM_SEED),
            input_fofn_fn,
            str(out_prefix)
            ]

    mod.main(argv)

    result_whitelist = open(str(out_prefix) + '.whitelist.json').read()
    result_all_zmw = open(str(out_prefix) + '.all.json').read()

    assert sorted(json.loads(exp_whitelist)) == sorted(json.loads(result_whitelist))
    assert sorted(json.loads(exp_all)) == sorted(json.loads(result_all_zmw))

@pytest.mark.parametrize('test_data', tests_strategy_random)
def test_run_strategy_random(tmpdir, test_data):
    check_run(tmpdir, **test_data)

@pytest.mark.parametrize('test_data', tests_strategy_longest)
def test_run_strategy_longest(tmpdir, test_data):
    check_run(tmpdir, **test_data)

@pytest.mark.parametrize('test_data', tests_strategy_random)
def test_main_random(tmpdir, test_data):
    check_main(tmpdir, **test_data)

@pytest.mark.parametrize('test_data', tests_strategy_longest)
def test_main_longest(tmpdir, test_data):
    check_main(tmpdir, **test_data)

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass
