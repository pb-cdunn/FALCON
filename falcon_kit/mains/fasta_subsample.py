#!/usr/bin/env python2.7
"""
Performs a single pass over an input FASTA/FOFN, and collects
all ZMWs. For each ZMW it calculates the expected molecular size by picking
the internal median subread length.
The script outputs a JSON file with a whitelist of ZMWs selected by a given
strategy (random, longest, etc.) and desired coverage of a genome.
Author: Ivan Sovic
"""
from falcon_kit.mains.fasta_filter import ZMWTuple
from falcon_kit.util.system import set_random_seed

import falcon_kit.FastaReader as FastaReader
import falcon_kit.mains.fasta_filter as fasta_filter
import falcon_kit.io as io

import os
import sys
import argparse
import logging
import contextlib
import itertools
import random
import json
import copy

LOG = logging.getLogger()

STRATEGY_RANDOM = 'random'
STRATEGY_LONGEST = 'longest'

def strategy_func_random(zmws):
    """
    >>> random.seed(12345); strategy_func_random([])
    []
    >>> random.seed(12345); strategy_func_random([('synthetic/1', 9)])
    [('synthetic/1', 9)]
    >>> random.seed(12345); strategy_func_random([('synthetic/1', 9), ('synthetic/2', 21), ('synthetic/3', 9), ('synthetic/4', 15), ('synthetic/5', 20)])
    [('synthetic/2', 21), ('synthetic/4', 15), ('synthetic/5', 20), ('synthetic/1', 9), ('synthetic/3', 9)]
    """
    ret = copy.deepcopy(zmws)
    random.shuffle(ret)
    return ret

def strategy_func_longest(zmws):
    """
    >>> strategy_func_longest([])
    []
    >>> strategy_func_longest([('synthetic/1', 9)])
    [('synthetic/1', 9)]
    >>> strategy_func_longest([('synthetic/1', 9), ('synthetic/2', 21), ('synthetic/3', 9), ('synthetic/4', 15), ('synthetic/5', 20)])
    [('synthetic/2', 21), ('synthetic/5', 20), ('synthetic/4', 15), ('synthetic/1', 9), ('synthetic/3', 9)]
    """
    return sorted(zmws, key = lambda x: x[1], reverse = True)

STRATEGY_TYPE_TO_FUNC = {   STRATEGY_RANDOM: strategy_func_random,
                            STRATEGY_LONGEST: strategy_func_longest
                        }

def get_strategy_func(strategy_type):
    """
    >>> get_strategy_func(STRATEGY_RANDOM) == strategy_func_random
    True
    >>> get_strategy_func(STRATEGY_LONGEST) == strategy_func_longest
    True
    >>> try:
    ...     get_strategy_func('nonexistent_strategy')
    ...     print('False')
    ... except:
    ...     print('True')
    True
    """
    assert strategy_type in STRATEGY_TYPE_TO_FUNC, 'Unknown strategy type: "{}"'.format(str(strategy_type))
    return STRATEGY_TYPE_TO_FUNC[strategy_type]

def select_zmws(zmws, min_requested_bases):
    """
    >>> select_zmws([], 0)
    ([], 0)
    >>> select_zmws([], 10)
    ([], 0)
    >>> select_zmws([('zmw/1', 1), ('zmw/2', 2), ('zmw/3', 5), ('zmw/4', 7), ('zmw/5', 10), ('zmw/6', 15)], 10)
    (['zmw/1', 'zmw/2', 'zmw/3', 'zmw/4'], 15)
    >>> select_zmws([('zmw/1', 1), ('zmw/2', 2), ('zmw/3', 5), ('zmw/4', 7), ('zmw/5', 10), ('zmw/6', 15)], 20)
    (['zmw/1', 'zmw/2', 'zmw/3', 'zmw/4', 'zmw/5'], 25)
    >>> select_zmws([('zmw/1', 1), ('zmw/1', 2), ('zmw/1', 5), ('zmw/1', 7), ('zmw/1', 10), ('zmw/1', 15)], 20)
    (['zmw/1', 'zmw/1', 'zmw/1', 'zmw/1', 'zmw/1'], 25)
    """
    # Select the first N ZMWs which sum up to the desired coverage.
    num_bases = 0
    subsampled_zmws = []
    for zmw_name, seq_len in zmws:
        num_bases += seq_len
        subsampled_zmws.append(zmw_name)
        if num_bases >= min_requested_bases:
            break
    return subsampled_zmws, num_bases

def calc_stats(total_unique_molecular_bases, total_bases, output_bases, genome_size, coverage):
    """
    >>> calc_stats(0, 0, 0, 0, 0) == \
    {'genome_size': 0, 'coverage': 0, 'total_bases': 0, 'total_unique_molecular_bases': 0, \
    'output_bases': 0, 'unique_molecular_avg_cov': 0.0, 'output_avg_cov': 0.0, 'total_avg_cov': 0.0}
    True
    >>> calc_stats(10000, 100000, 2000, 1000, 2) == \
    {'genome_size': 1000, 'coverage': 2, 'total_bases': 100000, 'total_unique_molecular_bases': 10000, \
    'output_bases': 2000, 'unique_molecular_avg_cov': 10.0, 'output_avg_cov': 2.0, 'total_avg_cov': 100.0}
    True
    """
    unique_molecular_avg_cov = 0.0 if genome_size == 0 else float(total_unique_molecular_bases) / float(genome_size)
    total_avg_cov = 0.0 if genome_size == 0 else float(total_bases) / float(genome_size)
    output_avg_cov = 0.0 if genome_size == 0 else float(output_bases) / float(genome_size)
    ret = {}
    ret['genome_size'] = genome_size
    ret['coverage'] = coverage
    ret['total_bases'] = total_bases
    ret['total_unique_molecular_bases'] = total_unique_molecular_bases
    ret['output_bases'] = output_bases
    ret['total_avg_cov'] = total_avg_cov
    ret['unique_molecular_avg_cov'] = unique_molecular_avg_cov
    ret['output_avg_cov'] = output_avg_cov
    return ret

def collect_zmws(yield_zmwtuple_func):
    """
    >>> collect_zmws([])
    ([], 0, 0)
    >>> collect_zmws([\
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=0, subread_end=1000, seq_len=1000, subread_record=None, subread_header='test/1/0_1000', subread_id=0), \
    ])
    ([('test/1', 1000)], 1000, 1000)
    >>> collect_zmws([\
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=123, subread_end=456, seq_len=1000, subread_record=None, subread_header='test/1/0_1000', subread_id=0), \
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=123, subread_end=456, seq_len=2000, subread_record=None, subread_header='test/1/1000_3000', subread_id=0), \
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=123, subread_end=456, seq_len=3000, subread_record=None, subread_header='test/1/3000_6000', subread_id=0), \
    ])
    ([('test/1', 2000)], 2000, 6000)
    >>> collect_zmws([\
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=123, subread_end=456, seq_len=1000, subread_record=None, subread_header='test/1/0_1000', subread_id=0), \
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=123, subread_end=456, seq_len=2000, subread_record=None, subread_header='test/1/1000_3000', subread_id=1), \
        ZMWTuple(movie_name='test' , zmw_id='1', subread_start=123, subread_end=456, seq_len=3000, subread_record=None, subread_header='test/1/3000_6000', subread_id=2), \
        ZMWTuple(movie_name='test' , zmw_id='2', subread_start=123, subread_end=456, seq_len=10000, subread_record=None, subread_header='header2', subread_id=3), \
    ])
    ([('test/1', 2000), ('test/2', 10000)], 12000, 16000)
    """
    zmws = []
    unique_molecular_size = 0
    total_size = 0
    for zmw_id, zmw_subreads in itertools.groupby(yield_zmwtuple_func, lambda x: x.zmw_id):
        zmw_subreads_list = list(zmw_subreads)
        zrec = fasta_filter.internal_median_zmw_subread(zmw_subreads_list)
        movie_zmw = zrec.movie_name + '/' + zrec.zmw_id
        unique_molecular_size += zrec.seq_len
        total_size += sum([zmw.seq_len for zmw in zmw_subreads_list])
        zmws.append((movie_zmw, zrec.seq_len))
    return zmws, unique_molecular_size, total_size

def yield_record(input_files):
    for input_fn in input_files:
        with open(input_fn, 'r') as fp_in:
            fasta_records = FastaReader.yield_fasta_record(fp_in, log=LOG.info)
            for record in fasta_records:
                yield record

def run(yield_zmw_tuple_func, coverage, genome_size, strategy_func):
    zmws, total_unique_molecular_bases, total_bases = collect_zmws(yield_zmw_tuple_func)
    zmws = strategy_func(zmws)
    subsampled_zmws, output_bases = select_zmws(zmws, coverage * genome_size)
    stats_dict = calc_stats(total_unique_molecular_bases, total_bases, output_bases, genome_size, coverage)
    return subsampled_zmws, zmws, stats_dict

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Produces a list of ZMW where the median unique molecular "\
                                        "coverage sums up to the desired coverage of the given genome size.s",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--strategy', type=str, default='random',
                        help='Subsampling strategy: random, longest')
    parser.add_argument('--coverage', type=float, default=60,
                        help='Desired coverage for subsampling.')
    parser.add_argument('--genome-size', type=float, default=0,
                        help='Genome size estimate of the input dataset.', required=True)
    parser.add_argument('--random-seed', type=int, default=12345,
                        help='Seed value used for the random generator.', required=False)
    parser.add_argument('input_fn', type=str, default='input.fofn',
                        help='Input FASTA files or a FOFN. (Streaming is not allowed).')
    parser.add_argument('out_prefix', type=str, default='input.cov',
                        help='Prefix of the output files to generate.')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)

    strategy_func = get_strategy_func(args.strategy)
    LOG.info('Using subsampling strategy: "{strategy}"'.format(strategy=args.strategy))

    set_random_seed(args.random_seed)

    input_files = list(io.yield_abspath_from_fofn(args.input_fn))

    zmws_whitelist, zmws_all, stats_dict  = run(
            fasta_filter.yield_zmwtuple(yield_record(input_files), None, False), args.coverage, args.genome_size, strategy_func)

    out_zmw_whitelist = args.out_prefix + '.whitelist.json'
    out_all_zmws = args.out_prefix + '.all.json'
    out_zmw_stats = args.out_prefix + '.stats.json'

    with open(out_zmw_whitelist, 'w') as fp_out_whitelist, \
         open(out_all_zmws, 'w') as fp_out_all_zmws, \
         open(out_zmw_stats, 'w') as fp_out_stats:

        # Write out the whitelist.
        fp_out_whitelist.write(json.dumps(zmws_whitelist))
        # Write the entire list of ZMWs and lengths, might be very informative.
        fp_out_all_zmws.write(json.dumps(zmws_all))
        # Write out the stats.
        fp_out_stats.write(json.dumps(stats_dict))

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
