from __future__ import absolute_import
from __future__ import print_function

import argparse
import sys
from falcon_kit import kup, falcon, DWA

class TooLongError(Exception): pass

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')

def get_aln_data(t_seq, q_seq):
    aln_data = []
    #x = []
    #y = []
    K = 8
    seq0 = t_seq
    lk_ptr = kup.allocate_kmer_lookup(1 << (K * 2))
    sa_ptr = kup.allocate_seq(len(seq0))
    sda_ptr = kup.allocate_seq_addr(len(seq0))
    kup.add_sequence(0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)
    q_id = "dummy"

    kmer_match_ptr = kup.find_kmer_pos_for_seq(
        q_seq, len(q_seq), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]

    if kmer_match.count != 0:
        aln_range_ptr = kup.find_best_aln_range(kmer_match_ptr, K, K * 5, 12)
        aln_range = aln_range_ptr[0]
        #x, y = list(zip(* [(kmer_match.query_pos[i], kmer_match.target_pos[i])
        #              for i in range(kmer_match.count)]))

        s1, e1, s2, e2 = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2

        log('Mapped (q, s1 = {}, e1 = {}, len1 = {}, (e1 - s1) = {}, t, s2 = {}, e2 = {}, (e2 - s2) = {}, len2 = {})'.format(
                s1, e1, e1 - s1, len(q_seq), s2, e2, e2 - s2, len(t_seq)))

        max_len = 250000 # to keep allocations < 16GB, given band_tol=1500
        if (e1 - s1) >= max_len or (e2 - s2) >= max_len:
            # DW.align() would crash, so raise here.
            # (500000 is the approx. upper bound for int overflow,
            #  but some users run out of memory anyway.)
            raise TooLongError('q_len={} or t_len={} are too big, over 500k'.format(
                (e1-s1), (e2-s2)))
        if e1 - s1 > 100:
            log('Calling DW_banded.align(q, s1 = {}, e1 = {}, len1 = {}, (e1 - s1) = {}, t, s2 = {}, e2 = {}, (e2 - s2) = {}, len2 = {})'.format(
                s1, e1, e1 - s1, len(q_seq), s2, e2, e2 - s2, len(t_seq)))
            alignment = DWA.align(q_seq[s1:e1], e1 - s1,
                                  seq0[s2:e2], e2 - s2,
                                  1500, 1)

            if alignment[0].aln_str_size > 100:
                aln_data.append((q_id, 0, s1, e1, len(q_seq), s2, e2, len(
                    seq0), alignment[0].aln_str_size, alignment[0].dist))
                aln_str1 = alignment[0].q_aln_str
                aln_str0 = alignment[0].t_aln_str

            DWA.free_alignment(alignment)

        kup.free_aln_range(aln_range_ptr)

    kup.free_kmer_match(kmer_match_ptr)
    kup.free_kmer_lookup(lk_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_seq_addr_array(sda_ptr)
    return aln_data #, x, y

def get_aln_results(ref_seq, query_seq, min_seq_len):
    # Align the a_ctg against the base.
    log('Aligning (DW): len(ref_seq) = %d, len(query_seq) = %d' % (len(ref_seq), len(query_seq)))
    delta_len = len(query_seq) - len(ref_seq)
    idt = 0.0
    cov = 0.0
    if len(ref_seq) > min_seq_len and len(query_seq) > min_seq_len:
        try:
            aln_data = get_aln_data(ref_seq, query_seq)
            if len(aln_data) != 0:
                idt = 1.0 - 1.0 * \
                    aln_data[-1][-1] / aln_data[-1][-2]
                cov = 1.0 * \
                    (aln_data[-1][3] - aln_data[-1]
                        [2]) / aln_data[-1][4]
            else:
                log('len(aln_data) == 0!')
        except TooLongError:
            log('WARNING: Seqs were too long for get_aln_data(), so we set idt/cov low enough to prevent filtering by dedup_a_tigs. len(ref_seq) = {}, len(query_seq) = {}'.format(len(ref_seq), len(query_seq)))
            idt = -1.0
            cov = -1.0
    return delta_len, idt, cov
