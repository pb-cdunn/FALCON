from __future__ import absolute_import
from __future__ import print_function

import argparse
import sys

from falcon_kit import kup, falcon
import edlib

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')

def count_cigar_ops(cigar):
    """
    For curious people: regexes are very slow for parsing CIGAR strings.
    """
    b = 0
    num_m, num_i, num_d = 0, 0, 0
    for i in xrange(len(cigar)):
        if cigar[i] <= '9':
            continue
        # Check if there are no digits before the op char.
        assert(b < i)
        count = int(cigar[b:i])
        op = cigar[i]
        b = i + 1
        if op == 'D':
            num_d += count
        elif op == 'I':
            num_i += count
        elif op in ['M', '=', 'X']:
            num_m += count
        else:        # pragma: no cover
            pass     # pragma: no cover
    # Check if there are dangling ops.
    assert(b == len(cigar))
    total_len = num_d + num_i + num_m
    return num_m, num_i, num_d, total_len

def get_aln_data(t_seq, q_seq):
    aln_data = []
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

        s1, e1, s2, e2 = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2

        log('Mapped (q, s1 = {}, e1 = {}, len1 = {}, (e1 - s1) = {}, t, s2 = {}, e2 = {}, (e2 - s2) = {}, len2 = {})'.format(
                s1, e1, e1 - s1, len(q_seq), s2, e2, e2 - s2, len(t_seq)))

        if e1 - s1 > 100:
            log('Calling edlib.align(q, s1 = {}, e1 = {}, len1 = {}, (e1 - s1) = {}, t, s2 = {}, e2 = {}, (e2 - s2) = {}, len2 = {})'.format(
                s1, e1, e1 - s1, len(q_seq), s2, e2, e2 - s2, len(t_seq)))

            # Align using Edlib instead of DWA.
            edlib_result = edlib.align(q_seq[s1:e1], seq0[s2:e2], mode="NW")

            delta_l = len(q_seq) - len(t_seq)
            cov = float(e1 - s1) / float(len(q_seq))
            idt = float(e1 - s1 - edlib_result['editDistance']) / float(e1 - s1)

            aln_data.append((q_id, 0, s1, e1, len(q_seq), s2, e2, len(seq0),
                            delta_l, idt, cov))

        kup.free_aln_range(aln_range_ptr)

    kup.free_kmer_match(kmer_match_ptr)
    kup.free_kmer_lookup(lk_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_seq_addr_array(sda_ptr)
    return aln_data #, x, y

def get_global_aln_results(ref_seq, query_seq, min_seq_len):
    """
    Aligns two sequences globally, and returns (delta_len, idt, cov) values
    compatible with the legacy deduplication code.
    Currently unused - it was used in an intermediate version, and might be useful
    at some point in the future.
    """
    log('Aligning (Edlib): len(ref_seq) = %d, len(query_seq) = %d' % (len(ref_seq), len(query_seq)))
    delta_len = len(query_seq) - len(ref_seq)

    idt = 0.0
    cov = 0.0

    if len(ref_seq) < min_seq_len or len(query_seq) < min_seq_len:
        return delta_len, idt, cov

    result = edlib.align(query_seq, ref_seq, mode="NW", task="path")

    cigar = result['cigar']
    num_m, num_i, num_d, total_len = count_cigar_ops(result['cigar'])
    num_eq = (num_m + num_i + num_d) - result['editDistance']
    num_x = num_m - num_eq

    idt_query = float(num_eq) / float(num_eq + num_x + num_i)
    idt_ref = float(num_eq) / float(num_eq + num_x + num_d)
    idt = min(idt_query, idt_ref)

    cov = 1.0

    log('  - Alignment stats: num_m = %d, num_i = %d, num_d = %d, total_len = %d, num_eq = %d, num_x = %d' % (num_m, num_i, num_d, total_len, num_eq, num_x))

    return delta_len, idt, cov

def get_aln_results(ref_seq, query_seq, min_seq_len):
    """
    Runs the legacy mapping code, and aligns the selected region using Edlib
    instead of the legacy DWA alignment with quadratic memory.
    """
    log('Aligning (Edlib): len(ref_seq) = %d, len(query_seq) = %d' % (len(ref_seq), len(query_seq)))

    delta_len = len(query_seq) - len(ref_seq)
    idt = 0.0
    cov = 0.0

    if len(ref_seq) < min_seq_len or len(query_seq) < min_seq_len:
        return delta_len, idt, cov

    aln_data = get_aln_data(ref_seq, query_seq)
    if len(aln_data) != 0:
        delta_len, idt, cov = aln_data[-1][8:11]

    else:
        log('len(aln_data) == 0!')

    return delta_len, idt, cov
