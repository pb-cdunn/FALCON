#! /usr/bin/env python2.7

import falcon_kit.FastaReader as FastaReader

import os
import sys
import argparse
import collections
import itertools
import logging

LOG = logging.getLogger()

ZMWTuple = collections.namedtuple('ZMWTuple', ['movie_name', 'zmw_id', 'subread_pos', 'seq_len', 'subread_record', 'subread_header', 'subread_id'])

def tokenize_header(seq_header):
    rid = seq_header.split()[0]
    movie_name, zmw_id, subread_pos = rid.split('/')
    return movie_name, zmw_id, subread_pos

def median_zmw_subread(zmw_subreads):
    """Return subread_record with median seq_len.
    zmw_subreads is a list of ZMWTuple.
    """
    assert len(zmw_subreads) != 0

    sorted_subreads = sorted(zmw_subreads, key = lambda x: x.seq_len)

    # Not really a median since we round to the floor, and the definition of
    # median value would be the average of the two middle elements.
    # However, we need a concrete subread associated with the median value.
    median_id = len(sorted_subreads) // 2

    return sorted_subreads[median_id]

##############################
### Streamed-median filter ###
##############################
def yield_zmwtuples(records):
    for record in records:
        movie_name, zmw_id, subread_pos = tokenize_header(record.name)
        zrec = ZMWTuple(movie_name=movie_name, zmw_id=zmw_id, subread_pos=subread_pos,
                        seq_len=len(record.sequence), subread_record=record,
                        subread_header=record.name, subread_id=0)
        yield zrec

def run_streamed_median(fp_in, fp_out, fn='-'):
    fasta_records = FastaReader.yield_fasta_records(fp_in, fn, log=LOG.info)
    for zmw_id, zmw_subreads in itertools.groupby(yield_zmwtuples(fasta_records), lambda x: x.zmw_id):
        median_zrec = median_zmw_subread(list(zmw_subreads))
        fp_out.write(str(median_zrec.subread_record))
        fp_out.write('\n')

##############################
### Pass filter.           ###
##############################
def run_pass_filter(fp_in, fp_out, fn):
    for record in FastaReader.yield_fasta_records(fp_in, fn, log=LOG.info):
        fp_out.write(str(record))
        fp_out.write('\n')

##################################
### Double-pass median filter. ###
##################################
def run_median_filter(fp_in, fp_out, fn):
    # Expect an actual file, not a stream.
    assert(os.path.exists(fn))

    # Needed to jump back for the second pass.
    fp_in_start = fp_in.tell()

    # Stores all subreads for a ZMW.
    zmw_dict = collections.defaultdict(list)

    # First pass, collect all ZMW info.
    for record in FastaReader.yield_fasta_records(fp_in, fn, log=LOG.info):
        movie_name, zmw_id, subread_pos = tokenize_header(record.name)
        # Store None instead of the actual record to free the memory after yield.
        zrec = ZMWTuple(movie_name=movie_name, zmw_id=zmw_id, subread_pos=subread_pos,
                        seq_len=len(record.sequence), subread_record=None,
                        subread_header=record.name, subread_id=len(zmw_dict[zmw_id]))
        zmw_dict[zmw_id].append(zrec)

    # For each ZMW, keep only one particular subread, specified by it's order of
    # appearance in the input FASTA file (stored in median_zrec.subread_id).
    whitelist = collections.defaultdict(int)
    for zmw_id, zmw_subreads in zmw_dict.iteritems():
        median_zrec = median_zmw_subread(list(zmw_subreads))
        whitelist[zmw_id] = median_zrec.subread_id

    # Second pass, yield selected sequences.
    # Rewind.
    fp_in.seek(fp_in_start, os.SEEK_SET)
    # Fly-through.
    for record in FastaReader.yield_fasta_records(fp_in, fn, log=LOG.info):
        movie_name, zmw_id, subread_pos = tokenize_header(record.name)
        # Write-out only one particular subread from the ZMW.
        if whitelist[zmw_id] == 0:
            fp_out.write(str(record))
            fp_out.write('\n')
        whitelist[zmw_id] -= 1

##############################
### Main and cmds.         ###
##############################
def cmd_run_pass_filter(args):
    if args.input_path == '-':
        fp_in = sys.stdin
    else:
        fp_in = open(args.input_path)

    run_pass_filter(fp_in, sys.stdout, args.input_path)

def cmd_run_streamed_median_filter(args):
    if args.input_path == '-':
        fp_in = sys.stdin
    else:
        fp_in = open(args.input_path)

    run_streamed_median(fp_in, sys.stdout, args.input_path)

def cmd_run_median_filter(args):
    with open(args.input_path, 'r') as fp_in:
        run_median_filter(fp_in, sys.stdout, args.input_path)

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_args(argv):
    description = 'Filters the input FASTA file according to one of the selected filters.'
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=HelpF,
    )

    subparsers = parser.add_subparsers(help='sub-command help')

    help_pass = 'The no-op filter - passes every FASTA record to stdout. If input_path is "-", input is read from stdin.'
    help_streamed_median = 'Applies the median-length ZMW filter by running a single-pass over the data. The input subreads should be groupped by ZMW. If input_path is "-", input is read from stdin.'
    help_median = 'Applies the median-length ZMW filter by running two passes over the data. Only one subread per ZMW is output, based on median-length selection. The input_path needs to be a file.'

    parser_pass = subparsers.add_parser('pass',
            formatter_class=HelpF,
            description=help_pass,
            help=help_pass)
    parser_pass.set_defaults(func=cmd_run_pass_filter)

    parser_median = subparsers.add_parser('median',
            formatter_class=HelpF,
            description=help_median,
            help=help_median)
    parser_median.set_defaults(func=cmd_run_median_filter)

    parser_median = subparsers.add_parser('streamed-median',
            formatter_class=HelpF,
            description=help_streamed_median,
            help=help_streamed_median)
    parser_median.set_defaults(func=cmd_run_streamed_median_filter)

    parser.add_argument('input_path', help='Input PacBio FASTA file')

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    args.func(args)

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
