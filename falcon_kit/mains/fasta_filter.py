#! /usr/bin/env python2.7

import falcon_kit.FastaReader as FastaReader

import os
import sys
import argparse
import collections
import itertools
import logging
import contextlib
import json

LOG = logging.getLogger()

ZMWTuple = collections.namedtuple('ZMWTuple', ['movie_name', 'zmw_id', 'subread_start', 'subread_end', 'seq_len', 'subread_record', 'subread_header', 'subread_id'])

def check_in_whitelist(whitelist_set, movie_name, zmw_id):
    """
    >>> check_in_whitelist([], 'foo', '1')
    False
    >>> check_in_whitelist(['bar/1'], 'foo', '1')
    False
    >>> check_in_whitelist(['foo/1'], 'foo', '1')
    True
    """
    movie_zmw = '{}/{}'.format(movie_name, zmw_id)
    return movie_zmw in whitelist_set

def tokenize_header(seq_header):
    """
    >>> tokenize_header('foo/123/0_100')
    ('foo', '123', 0, 100)
    """
    rid = seq_header.split()[0]
    movie_name, zmw_id, subread_pos = rid.split('/')
    subread_start, subread_end = [int(val) for val in subread_pos.split('_')]
    return movie_name, zmw_id, subread_start, subread_end

def yield_record_and_tokenized_headers(whitelist_set, records):
    """For each record, yield (record, tokens)
    but only if whitelisted.

    records: iterable
    whitelist_set: has __contains__, but empty means use everything.
    """
    for record in records:
        tokens = tokenize_header(record.name)
        movie_name, zmw_id, subread_start, subread_end = tokens
        if whitelist_set and not check_in_whitelist(whitelist_set, movie_name, zmw_id):
            continue
        yield record, tokens

def yield_record(whitelist_set, records):
    """Yield each record,
    but only if whitelisted.

    This is an optimized version of yield_record_and_tokenized_headers(),
    to avoid tokenizing when we have no whitelist.

    records: iterable
    whitelist_set: has __contains__, but empty means use everything.
    """
    for record in records:
        if not whitelist_set:
            # no need to tokenize
            yield record
            continue
        tokens = tokenize_header(record.name)
        movie_name, zmw_id, subread_start, subread_end = tokens
        if not check_in_whitelist(whitelist_set, movie_name, zmw_id):
            continue
        yield record

def longest_zmw_subread(zmw_subreads):
    """Return subread_record with longest seq_len.
    zmw_subreads is a list of ZMWTuple.
    """
    assert len(zmw_subreads) != 0

    return max(zmw_subreads, key = lambda x: x.seq_len)

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

def internal_median_zmw_subread(zmw_subreads):
    """Returns a single subread based on the following selection criteria:
    - If the ZMW has < 3 subreads, the maximum one is output.
    - If the ZMW has >= 3 subreads, the median one is selected only from the internal
      ones (ignoring the first and the last subread)
    This is intended to prevent the impact of very short first and last subread on the
    median selection, since the polymerase can start/stop in the middle of the insert.
    zmw_subreads is a list of ZMWTuple.
    """
    assert len(zmw_subreads) != 0

    selected_subread = None

    if len(zmw_subreads) < 3:
        sorted_subreads = sorted(zmw_subreads, key = lambda x: x.seq_len)
        selected_subread = sorted_subreads[-1]
    else:
        sorted_by_pos = sorted(zmw_subreads, key = lambda x: x.subread_start)
        sorted_subreads = sorted(sorted_by_pos[1:-1], key = lambda x: x.seq_len)
        median_id = len(sorted_subreads) // 2
        selected_subread = sorted_subreads[median_id]

    return selected_subread

##############################
### Streamed-median filter ###
##############################
def yield_zmwtuple(records, whitelist_set, store_record):
    subread_id = 0
    for (record, tokens) in yield_record_and_tokenized_headers(whitelist_set, records):
        movie_name, zmw_id, subread_start, subread_end = tokens
        record_to_store = record if store_record else None
        zrec = ZMWTuple(movie_name=movie_name, zmw_id=zmw_id,
                        subread_start=subread_start, subread_end=subread_end,
                        seq_len=len(record.sequence), subread_record=record_to_store,
                        subread_header=record.name, subread_id=subread_id)
        subread_id += 1
        yield zrec

def write_streamed(fp_out, yield_zmwtuple_func, zmw_filter_func):
    for zmw_id, zmw_subreads in itertools.groupby(yield_zmwtuple_func(store_record=True), lambda x: x.zmw_id):
        zrec = zmw_filter_func(list(zmw_subreads))
        fp_out.write(str(zrec.subread_record))
        fp_out.write('\n')

def run_streamed_median_filter(fp_in, fp_out, whitelist_set, zmw_filter_func=median_zmw_subread):
    def yield_zmwtuple_func(store_record=True):
        fasta_records = FastaReader.yield_fasta_record(fp_in, log=LOG.info)
        return yield_zmwtuple(fasta_records, whitelist_set, store_record)
    write_streamed(fp_out, yield_zmwtuple_func, zmw_filter_func)

##############################
### Longest filter.
##############################
def run_streamed_longest_filter(fp_in, fp_out, whitelist_set):
    def yield_zmwtuple_func(store_record=True):
        fasta_records = FastaReader.yield_fasta_record(fp_in, log=LOG.info)
        return yield_zmwtuple(fasta_records, whitelist_set, store_record)
    write_streamed(fp_out, yield_zmwtuple_func, zmw_filter_func=longest_zmw_subread)

##############################
### Pass filter.           ###
##############################
def run_pass_filter(fp_in, fp_out, whitelist_set):
    for record in yield_record(whitelist_set, FastaReader.yield_fasta_record(fp_in, log=LOG.info)):
        fp_out.write(str(record))
        fp_out.write('\n')

##################################
### Double-pass median filter. ###
##################################
def write_doublepass_median(fp_out, yield_zmwtuple_func, zmw_filter_func=median_zmw_subread):
    # Stores all subreads for a ZMW.
    zmw_dict = collections.defaultdict(list)

    # First pass, collect all ZMW info.
    for zrec in yield_zmwtuple_func(store_record=False):
        # Store None instead of the actual record to free the memory after yield.
        zmw_id = zrec.zmw_id
        zmw_dict[zmw_id].append(zrec)

    # For each ZMW, keep only one particular subread.
    selected = collections.defaultdict(int)
    for zmw_id, zmw_subreads in zmw_dict.iteritems():
        median_zrec = zmw_filter_func(list(zmw_subreads))
        selected[zmw_id] = median_zrec.subread_id

    # Second pass, yield selected sequences.
    # This must be exactly the same order, so we end up with the same computed subread_id for each record.
    for zrec in yield_zmwtuple_func(store_record=True):
        zmw_id = zrec.zmw_id
        subread_id = zrec.subread_id
        if selected[zmw_id] == subread_id:
            fp_out.write(str(zrec.subread_record))
            fp_out.write('\n')

def run_median_filter(fp_in, fp_out, whitelist_set, zmw_filter_func=median_zmw_subread):
    # Needed to jump back for the second pass.
    # Expect an actual file, not a pipe.
    try:
        fp_in_start = fp_in.tell()
        fp_in.seek(fp_in_start)
    except Exception as e:
        msg = 'fileobj.tell()/seek() failed. Cannot rewind, so cannot do multi-pass. {}\n{}'.format(
            type(fp_in), dir(fp_in))
        raise AssertionError(msg, e)

    def yield_zmwtuple_func(store_record=True):
        # Rewind.
        fp_in.seek(fp_in_start, os.SEEK_SET)

        fasta_records = FastaReader.yield_fasta_record(fp_in, log=LOG.info)
        return yield_zmwtuple(fasta_records, whitelist_set, store_record)

    write_doublepass_median(fp_out, yield_zmwtuple_func, zmw_filter_func)

###############################
### Internal median filter. ###
###############################
def run_internal_median_filter(fp_in, fp_out, whitelist_set):
    run_median_filter(fp_in, fp_out, whitelist_set, zmw_filter_func=internal_median_zmw_subread)

#######################################
### Streamed internal median filter ###
#######################################
def run_streamed_internal_median_filter(fp_in, fp_out, whitelist_set):
    run_streamed_median_filter(fp_in, fp_out, whitelist_set=whitelist_set, zmw_filter_func=internal_median_zmw_subread)

##############################
### Main and cmds.         ###
##############################
@contextlib.contextmanager
def open_input_stream(input_path):
    """stdin if '-'
    """
    if input_path == '-':
        yield sys.stdin
    else:
        with open(input_path) as stream:
            yield stream

def load_zmw_whitelist(zmw_whitelist_fn):
    """Read from json filename, or do nothing if null fn.
    Return as a set, empty-set if empty.
    Raise if missing non-null fn.
    """
    ret = set()
    if not zmw_whitelist_fn:
        return set()
    with open(zmw_whitelist_fn, 'r') as fp_in:
        try:
            return set(json.loads(fp_in.read()))
        except ValueError:
            LOG.error('Failed to parse "{}" as JSON. Assuming empty whitelist.'.format(zmw_whitelist_fn))
            return set()

def cmd_run_pass_filter(args):
    whitelist_set = load_zmw_whitelist(args.zmw_whitelist_fn)
    with open_input_stream(args.input_path) as fp_in:
        run_pass_filter(fp_in, sys.stdout, whitelist_set)

def cmd_run_streamed_longest_filter(args):
    whitelist_set = load_zmw_whitelist(args.zmw_whitelist_fn)
    with open_input_stream(args.input_path) as fp_in:
        run_streamed_longest_filter(fp_in, sys.stdout, whitelist_set)

def cmd_run_streamed_median_filter(args):
    whitelist_set = load_zmw_whitelist(args.zmw_whitelist_fn)
    with open_input_stream(args.input_path) as fp_in:
        run_streamed_median_filter(fp_in, sys.stdout, whitelist_set)

def cmd_run_median_filter(args):
    whitelist_set = load_zmw_whitelist(args.zmw_whitelist_fn)
    # Don't allow '-' for the double-pass median filter.
    with open(args.input_path, 'r') as fp_in:
        run_median_filter(fp_in, sys.stdout, whitelist_set)

def cmd_run_internal_median_filter(args):
    whitelist_set = load_zmw_whitelist(args.zmw_whitelist_fn)
    with open(args.input_path, 'r') as fp_in:
        run_internal_median_filter(fp_in, sys.stdout, whitelist_set)

def cmd_run_streamed_internal_median_filter(args):
    whitelist_set = load_zmw_whitelist(args.zmw_whitelist_fn)
    with open_input_stream(args.input_path) as fp_in:
        run_streamed_internal_median_filter(fp_in, sys.stdout, whitelist_set)

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_args(argv):
    description = 'Filters the input FASTA file according to one of the selected filters.'
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=HelpF,
    )

    parser.add_argument(
        '--zmw-whitelist-fn', default='',
        help='A JSON file containing a list of "movie_name/zmw" IDs to retain. If the file is an empty list, all ZMWs will be used; otherwise, only the listed ones will be whitelisted. (Applies to all filters.)',
        required = False,
    )

    subparsers = parser.add_subparsers(help='sub-command help')

    help_pass = 'The no-op filter - passes every FASTA record to stdout. If input_path is "-", input is read from stdin.'
    help_streamed_longest = 'Selects the longest read in each ZMW by running a single-pass over the data (i.e. "streamed"). The input subreads should be groupped by ZMW. If input_path is "-", input is read from stdin.'
    help_median = 'Applies the median-length ZMW filter by running two passes over the data. Only one subread per ZMW is output, based on median-length selection. The input_path needs to be a file.'
    help_streamed_median = 'Applies the median-length ZMW filter by running a single-pass over the data. The input subreads should be groupped by ZMW. If input_path is "-", input is read from stdin.'
    help_internal_median = 'Applies the median-length ZMW filter only on internal subreads (ZMWs with >= 3 subreads) by running two passes over the data. For ZMWs with < 3 subreads, the maximum-length one is selected. The input_path needs to be a file.'
    help_streamed_internal_median = 'Applies the median-length ZMW filter only on internal subreads (ZMWs with >= 3 subreads) by running a single pass over the data. The input subreads should be groupped by ZMW. For ZMWs with < 3 subreads, the maximum-length one is selected. If input_path is "-", input is read from stdin.'

    parser_pass = subparsers.add_parser('pass',
            formatter_class=HelpF,
            description=help_pass,
            help=help_pass)
    parser_pass.set_defaults(func=cmd_run_pass_filter)

    parser_streamed_longest = subparsers.add_parser('streamed-longest',
            formatter_class=HelpF,
            description=help_streamed_longest,
            help=help_streamed_longest)
    parser_streamed_longest.set_defaults(func=cmd_run_streamed_longest_filter)

    parser_median = subparsers.add_parser('median',
            formatter_class=HelpF,
            description=help_median,
            help=help_median)
    parser_median.set_defaults(func=cmd_run_median_filter)

    parser_streamed_median = subparsers.add_parser('streamed-median',
            formatter_class=HelpF,
            description=help_streamed_median,
            help=help_streamed_median)
    parser_streamed_median.set_defaults(func=cmd_run_streamed_median_filter)

    parser_internal_median = subparsers.add_parser('internal-median',
            formatter_class=HelpF,
            description=help_internal_median,
            help=help_internal_median)
    parser_internal_median.set_defaults(func=cmd_run_internal_median_filter)

    parser_streamed_internal_median = subparsers.add_parser('streamed-internal-median',
            formatter_class=HelpF,
            description=help_streamed_internal_median,
            help=help_streamed_internal_median)
    parser_streamed_internal_median.set_defaults(func=cmd_run_streamed_internal_median_filter)

    parser.add_argument('input_path', help='Input PacBio FASTA file')

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    args.func(args)

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
