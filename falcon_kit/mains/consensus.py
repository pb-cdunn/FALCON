from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from builtins import range
from ctypes import (POINTER, c_char_p, c_uint, c_uint,
                    c_uint, c_uint, c_uint, c_double, string_at, pointer)
from falcon_kit.multiproc import Pool
from falcon_kit import falcon
import argparse
import logging
import multiprocessing
import os
import re
import sys
import falcon_kit
import falcon_kit.util.io as io
import collections

LOG = logging.getLogger()

falcon.generate_consensus.argtypes = [
    POINTER(c_char_p), c_uint, c_uint, c_uint, c_double]
falcon.generate_consensus.restype = POINTER(falcon_kit.ConsensusData)
falcon.free_consensus_data.argtypes = [POINTER(falcon_kit.ConsensusData)]

falcon.generate_consensus_from_mapping.argtypes = [
    POINTER(c_char_p), POINTER(POINTER(falcon_kit.AlnRange)), c_uint, c_uint, c_uint, c_double]
falcon.generate_consensus_from_mapping.restype = POINTER(falcon_kit.ConsensusData)

"""
SeqTuple encodes a single line in a block for consensus. Legacy code used only the 'name' and 'seq' (read from input),
but if the coordinates are already known, we can use this info.
The `qlen` and `tlen` are necessary because this consensus code can clip the end of a sequence if it's
beyond a certain threshold. If it's clipped, the start/end coordinates can fall within the clipped region,
which means that the internal alignment will have to be triggered.
The 'tstart' and 'tend' relate to the seed read, and the 'qstart' and 'qend' to the currenr read on the same line.
The current query should be in the same strand as the target. For consistency, we added a 'qstranq' as well, but
in the current LA4Falcon output it will always be 0.
The 'aln' field can be used to provide an alignment directly from the tool which determined that these sequences
need to go into the same block. This could be used downstream to prevent quadratic memory consumption during error
correction, and speed up the process.
Parameter 'is_trimmed' is a bool, indicating that the sequence was trimmed from the back because it exceeded the maximum length.
"""
SeqTuple = collections.namedtuple('SeqTuple', ['name', 'seq', 'qstrand', 'qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen', 'aln', 'is_mapped', 'is_trimmed'])

def get_longest_reads(seqs, max_n_read, max_cov_aln, sort=True):
    # including the sort kwarg allows us to avoid a redundant sort
    # in get_consensus_trimmed()
    if sort:
        seqs = seqs[:1] + sorted(seqs[1:], key=lambda x: -len(x.seq))

    longest_n_reads = max_n_read
    if max_cov_aln > 0:
        longest_n_reads = 1
        seed_len = len(seqs[0].seq)
        read_cov = 0
        for seq in seqs[1:]:
            if read_cov // seed_len > max_cov_aln:
                break
            longest_n_reads += 1
            read_cov += len(seq.seq)

        longest_n_reads = min(longest_n_reads, max_n_read)

    return(seqs[:longest_n_reads])


def get_alignment(seq1, seq0, edge_tolerance=1000):

    kup = falcon_kit.kup
    K = 8
    lk_ptr = kup.allocate_kmer_lookup(1 << (K * 2))
    sa_ptr = kup.allocate_seq(len(seq0))
    sda_ptr = kup.allocate_seq_addr(len(seq0))
    kup.add_sequence(0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)

    kup.mask_k_mer(1 << (K * 2), lk_ptr, 16)
    kmer_match_ptr = kup.find_kmer_pos_for_seq(
        seq1, len(seq1), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]
    aln_range_ptr = kup.find_best_aln_range2(kmer_match_ptr, K, K * 50, 25)
    #x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count )] )
    aln_range = aln_range_ptr[0]
    kup.free_kmer_match(kmer_match_ptr)
    s1, e1, s0, e0, km_score = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2, aln_range.score
    e1 += K + K // 2
    e0 += K + K // 2
    kup.free_aln_range(aln_range)
    len_1 = len(seq1)
    len_0 = len(seq0)
    if e1 > len_1:
        e1 = len_1
    if e0 > len_0:
        e0 = len_0

    aln_size = 1
    if e1 - s1 > 500:

        aln_size = max(e1 - s1, e0 - s0)
        aln_score = int(km_score * 48)
        aln_q_s = s1
        aln_q_e = e1
        aln_t_s = s0
        aln_t_e = e0

    kup.free_seq_addr_array(sda_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_kmer_lookup(lk_ptr)

    if s1 > edge_tolerance and s0 > edge_tolerance:
        return 0, 0, 0, 0, 0, 0, "none"

    if len_1 - e1 > edge_tolerance and len_0 - e0 > edge_tolerance:
        return 0, 0, 0, 0, 0, 0, "none"

    if e1 - s1 > 500 and aln_size > 500:
        return s1, s1 + aln_q_e - aln_q_s, s0, s0 + aln_t_e - aln_t_s, aln_size, aln_score, "aln"
    else:
        return 0, 0, 0, 0, 0, 0, "none"


def get_trimmed_seq(seq, s, e):
    # Mapping info is useless after clipping, so just reset it.
    ret = SeqTuple(name = seq.name, seq = seq.seq[s:e],
                    qstrand = seq.qstrand, qstart = -1, qend = -1, qlen = -1,
                    tstart = -1, tend = -1, tlen = -1,
                    aln = '*', is_mapped = False, is_trimmed = True)
    return ret

def get_consensus_core(seqs, min_cov, K, min_idt, allow_external_mapping):
    seqs_ptr = (c_char_p * len(seqs))()
    seqs_ptr[:] = [val.seq for val in seqs]

    all_seqs_mapped = False

    if allow_external_mapping:
        all_seqs_mapped = True
        for seq in seqs:
            if not seq.is_mapped:
                all_seqs_mapped = False
                break

    if not all_seqs_mapped:
        LOG.info('Internally mapping the sequences.')
        consensus_data_ptr = falcon.generate_consensus(
            seqs_ptr, len(seqs), min_cov, K, min_idt)

    else:
        LOG.info('Using external mapping coordinates from input.')
        aln_ranges_ptr = (POINTER(falcon_kit.AlnRange) * len(seqs))()
        for i, seq in enumerate(seqs):
            a = falcon_kit.AlnRange(seq.qstart, seq.qend, seq.tstart, seq.tend, (seq.qend - seq.qstart))
            aln_ranges_ptr[i] = pointer(a)
        consensus_data_ptr = falcon.generate_consensus_from_mapping(
            seqs_ptr, aln_ranges_ptr, len(seqs), min_cov, K, min_idt)
        del aln_ranges_ptr

    del seqs_ptr

    if not consensus_data_ptr:
        return ''
    # assert consensus_data_ptr
    consensus = string_at(consensus_data_ptr[0].sequence)[:]
    eff_cov = consensus_data_ptr[0].eff_cov[:len(consensus)]
    LOG.debug(' Freeing')
    falcon.free_consensus_data(consensus_data_ptr)
    return consensus

def get_consensus_without_trim(c_input):
    seqs, seed_id, config = c_input
    LOG.debug('Starting get_consensus_without_trim(len(seqs)=={}, seed_id={})'.format(
        len(seqs), seed_id))
    min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln, allow_external_mapping = config
    if len(seqs) > max_n_read:
        seqs = get_longest_reads(seqs, max_n_read, max_cov_aln, sort=True)

    consensus = get_consensus_core(seqs, min_cov, K, min_idt, allow_external_mapping)
    LOG.debug(' Finishing get_consensus_without_trim(seed_id={})'.format(seed_id))

    return consensus, seed_id

def get_consensus_with_trim(c_input):
    seqs, seed_id, config = c_input
    LOG.debug('Starting get_consensus_with_trim(len(seqs)=={}, seed_id={})'.format(
        len(seqs), seed_id))
    min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln, allow_external_mapping = config
    trim_seqs = []
    seed = seqs[0]
    for seq in seqs[1:]:
        aln_data = get_alignment(seq.seq, seed.seq, edge_tolerance)
        s1, e1, s2, e2, aln_size, aln_score, c_status = aln_data
        if c_status == "none":
            continue
        if aln_score > 1000 and e1 - s1 > 500:
            e1 -= trim_size
            s1 += trim_size
            trim_seqs.append((e1 - s1, get_trimmed_seq(seq, s1, e1)))
            # trim_seqs.append((e1 - s1, seq.seq[s1:e1]))
    trim_seqs.sort(key=lambda x: -x[0])  # use longest alignment first
    trim_seqs = [x[1] for x in trim_seqs]

    trim_seqs = [seed] + trim_seqs
    if len(trim_seqs[1:]) > max_n_read:
        # seqs already sorted, dont' sort again
        trim_seqs = get_longest_reads(
            trim_seqs, max_n_read, max_cov_aln, sort=False)

    consensus = get_consensus_core(trim_seqs, min_cov, K, min_idt, allow_external_mapping)
    LOG.debug(' Finishing get_consensus_with_trim(seed_id={})'.format(seed_id))

    return consensus, seed_id

def get_seq_data(config, min_n_read, min_len_aln):
    max_len = 128000
    min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln, allow_external_mapping = config
    seqs = []
    seed_id = None
    seed_len = 0
    seqs_data = []
    read_cov = 0
    read_ids = set()
    with sys.stdin as f:
        for line in f:
            split_line = line.strip().split()
            if len(split_line) < 2:
                continue

            qname = split_line[0]
            qseq = split_line[1]
            qstrand, qstart, qend, qlen = 0, -1, -1, -1
            tstart, tend, tlen = -1, -1, -1
            aln, is_mapped, is_trimmed = '*', False, False

            if len(split_line) >= 10:
                qstrand = int(split_line[2])
                qstart = int(split_line[3])
                qend = int(split_line[4])
                qlen = int(split_line[5])
                tstart = int(split_line[6])
                tend = int(split_line[7])
                tlen = int(split_line[8])
                aln = split_line[9]
                is_mapped = True

            new_seq = SeqTuple(name = qname, seq = qseq,
                                qstrand = qstrand, qstart = qstart, qend = qend, qlen = qlen,
                                tstart = tstart, tend = tend, tlen = tlen,
                                aln = aln, is_mapped = is_mapped, is_trimmed = is_trimmed)

            if len(new_seq.seq) > max_len:
                new_seq = get_trimmed_seq(new_seq, 0, max_len - 1)

            if new_seq.name not in ("+", "-", "*"):
                if len(new_seq.seq) >= min_len_aln:
                    if len(seqs) == 0:
                        seqs.append(new_seq)  # the "seed"
                        seed_len = len(new_seq.seq)
                        seed_id = new_seq.name
                    if new_seq.name not in read_ids:  # avoidng using the same read twice. seed is used again here by design
                        seqs.append(new_seq)
                        read_ids.add(new_seq.name)
                        read_cov += len(new_seq.seq)

            elif split_line[0] == "+":
                if len(seqs) >= min_n_read and read_cov // seed_len >= min_cov_aln:
                    seqs = get_longest_reads(
                        seqs, max_n_read, max_cov_aln, sort=True)
                    yield (seqs, seed_id, config)
                #seqs_data.append( (seqs, seed_id) )
                seqs = []
                read_ids = set()
                seed_id = None
                read_cov = 0
            elif split_line[0] == "*":
                seqs = []
                read_ids = set()
                seed_id = None
                read_cov = 0
            elif split_line[0] == "-":
                # yield (seqs, seed_id)
                #seqs_data.append( (seqs, seed_id) )
                break


def format_seq(seq, col):
    return "\n".join([seq[i:(i + col)] for i in range(0, len(seq), col)])


def parse_args(argv):
    parser = argparse.ArgumentParser(description='a simple multi-processor consensus sequence generator',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n-core', type=int, default=24,
                        help='number of processes used for generating consensus; '
                        '0 for main process only')
    parser.add_argument('--min-cov', type=int, default=6,
                        help='minimum coverage to break the consensus')
    parser.add_argument('--min-cov-aln', type=int, default=10,
                        help='minimum coverage of alignment data; a seed read with less than MIN_COV_ALN average depth' +
                        ' of coverage will be completely ignored')
    parser.add_argument('--max-cov-aln', type=int, default=0,  # 0 to emulate previous behavior
                        help='maximum coverage of alignment data; a seed read with more than MAX_COV_ALN average depth' + \
                        ' of coverage of the longest alignments will be capped, excess shorter alignments will be ignored')
    parser.add_argument('--min-len-aln', type=int, default=0,  # 0 to emulate previous behavior
                        help='minimum length of a sequence in an alignment to be used in consensus; any shorter sequence will be completely ignored')
    parser.add_argument('--min-n-read', type=int, default=10,
                        help='1 + minimum number of reads used in generating the consensus; a seed read with fewer alignments will ' +
                        'be completely ignored')
    parser.add_argument('--max-n-read', type=int, default=500,
                        help='1 + maximum number of reads used in generating the consensus')
    parser.add_argument('--trim', action="store_true", default=False,
                        help='trim the input sequence with k-mer spare dynamic programming to find the mapped range')
    parser.add_argument('--output-full', action="store_true", default=False,
                        help='output uncorrected regions too')
    parser.add_argument('--output-multi', action="store_true", default=False,
                        help='output multi correct regions')
    parser.add_argument('--min-idt', type=float, default=0.70,
                        help='minimum identity of the alignments used for correction')
    parser.add_argument('--edge-tolerance', type=int, default=1000,
                        help='for trimming, the there is unaligned edge leng > edge_tolerance, ignore the read')
    parser.add_argument('--trim-size', type=int, default=50,
                        help='the size for triming both ends from initial sparse aligned region')
    parser.add_argument('--allow-external-mapping', action="store_true", default=False,
                        help='if provided, externally determined mapping coordinates will be used for error correction')
    parser.add_argument('-v', '--verbose-level', type=float, default=2.0,
                        help='logging level (WARNING=3, INFO=2, DEBUG=1)')
    return parser.parse_args(argv[1:])

def run(args):
    logging.basicConfig(level=int(round(10*args.verbose_level)))

    assert args.n_core <= multiprocessing.cpu_count(), 'Requested n_core={} > cpu_count={}'.format(
            args.n_core, multiprocessing.cpu_count())

    def Start():
        LOG.info('Started a worker in {} from parent {}'.format(
            os.getpid(), os.getppid()))
    exe_pool = Pool(args.n_core, initializer=Start)
    if args.trim:
        get_consensus = get_consensus_with_trim
    else:
        get_consensus = get_consensus_without_trim

    K = 8
    config = args.min_cov, K, \
        args.max_n_read, args.min_idt, args.edge_tolerance, \
        args.trim_size, args.min_cov_aln, args.max_cov_aln, \
        args.allow_external_mapping
    # TODO: pass config object, not tuple, so we can add fields
    inputs = []
    for datum in get_seq_data(config, args.min_n_read, args.min_len_aln):
        inputs.append((get_consensus, datum))
    try:
        LOG.info('running {!r}'.format(get_consensus))
        for res in exe_pool.imap(io.run_func, inputs):
            process_get_consensus_result(res, args)
        LOG.info('finished {!r}'.format(get_consensus))
    except:
        LOG.exception('failed gen_consensus')
        exe_pool.terminate()
        raise

good_region = re.compile("[ACGT]+")

def process_get_consensus_result(res, args, limit=500):
        cns, seed_id = res
        seed_id = int(seed_id)
        if len(cns) < limit:
            return

        if args.output_full:
            print('>{:d}_f'.format(seed_id))
            print(cns)
        else:
            cns = good_region.findall(cns)
            if args.output_multi:
                seq_i = 0
                for cns_seq in cns:
                    if len(cns_seq) < limit:
                        continue
                    if seq_i >= 10:
                        break
                    print(">prolog/%s%01d/%d_%d" % (seed_id, seq_i, 0, len(cns_seq)))
                    print(format_seq(cns_seq, 80))
                    seq_i += 1
            else:
                if len(cns) == 0:
                    return
                cns.sort(key=lambda x: len(x))
                print('>{:d}'.format(seed_id))
                print(cns[-1])

def main(argv=sys.argv):
    args = parse_args(argv)
    run(args)

if __name__ == "__main__":
    main(sys.argv)
