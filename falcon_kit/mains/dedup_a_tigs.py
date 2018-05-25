from __future__ import absolute_import
from __future__ import print_function

from falcon_kit.FastaReader import open_fasta_reader
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

        max_len = 250000 # to keep allocations < 16GB, given band_tol=1500
        if (e1 - s1) >= max_len or (e2 - s2) >= max_len:
            # DW.align() would crash, so raise here.
            # (500000 is the approx. upper bound for int overflow,
            #  but some users run out of memory anyway.)
            raise TooLongError('q_len={} or t_len={} are too big, over 500k'.format(
                (e1-s1), (e2-s2)))
        if e1 - s1 > 100:
            log('Calling DW_banded.align(q, {}, t, {}, 1500, 1)'.format(
                e1-s1, e2-s2))
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

def get_aln_results(ref_seq, query_seq):
    # Align the a_ctg against the base.
    log('Aligning: len(ref_seq) = %d, len(query_seq) = %d' % (len(ref_seq), len(query_seq)))
    delta_len = len(query_seq) - len(ref_seq)
    idt = 0.0
    cov = 0.0
    if len(ref_seq) > 2000 and len(query_seq) > 2000:
        try:
            aln_data = get_aln_data(ref_seq, query_seq)
            if len(aln_data) != 0:
                idt = 1.0 - 1.0 * \
                    aln_data[-1][-1] / aln_data[-1][-2]
                cov = 1.0 * \
                    (aln_data[-1][3] - aln_data[-1]
                        [2]) / aln_data[-1][4]
        except TooLongError:
            log('WARNING: Seqs were too long for get_aln_data(), so we set idt/cov low enough to prevent filtering by dedup_a_tigs. len(ref_seq) = %d, len(query_seq) = %d'.format(len(ref_seq), len(query_seq)))
            idt = -1.0
            cov = -1.0
    return delta_len, idt, cov

def yield_single_compound(fp_in):
    """
    Loads all a_ctg with the same ctg_id and a_id.
    The header of a_ctg is:
        p_ctg_id-a_id-sub_id v w total_length total_score num_path_edges delta_len idt cov
    The first sequence in the returned list is the base a_ctg (the exact sequence
    which is part of the primary contig (the base of the bubble)).
    """
    ret = []
    prev_id = None
    for r in fp_in:
        tig_id, v, w, len_, ovl, ne, delta_l, idt, cov = r.name.split()
        p_ctg_id, a_id, sub_id = tig_id.split('-')
        curr_id = (p_ctg_id, a_id)
        if prev_id != None and prev_id != curr_id:
            yield ret
            ret = []
        prev_id = curr_id
        ret.append(r)
    yield ret

def filter_duplicate(compound_a_ctg, max_idt, max_aln_cov, min_len_diff, ploidy):
    """
    Takes a list of a_ctg sequences in a compound unitig (bubble) which need
    to be deduplicated according to the parameters.
    The zeroth sequence in the list is the "base" sequence. This sequence is
    already part of the primary path by definition, and will not be output.
    """

    # Sanity check.
    if len(compound_a_ctg) == 0:
        return []

    ret = []
    ref_seqs = [compound_a_ctg[0]]

    # Zeroth sequence is the base seq.
    for i in xrange(1, len(compound_a_ctg)):
        header = compound_a_ctg[i].name.split()
        a_ctg_id, v, w, len_, ovl, ne, delta_l, idt, cov = header
        a_ctg_seq = compound_a_ctg[i].sequence

        # Reset the values
        delta_l, idt, cov = 0.0, 1.0, 1.0
        is_duplicate = False

        # Align against the base sequence and all non-filtered alternate branches.
        loop_to = len(ref_seqs) if ploidy <= 0 else min(ploidy, len(ref_seqs))
        for j in xrange(0, loop_to):
            # Just fetch the components for readibility.
            ref_ctg_id = ref_seqs[j].name.split()[0]
            ref_seq = ref_seqs[j].sequence

            log('[i = %d, j = %d] Comparing: query "%s" vs ref "%s".' % (i, j, a_ctg_id, ref_ctg_id))

            # Align.
            delta_l, idt, cov = get_aln_results(ref_seq, a_ctg_seq)

            # Round to the floor of 2 decimal places. Needed to reproduce
            # old behaviour.
            idt = float('%.2f' % (idt))
            cov = float('%.2f' % (cov))

            log('  Rounded: new_delta_l = %d, new_idt = %.2f, new_cov = %.2f' % (delta_l, idt, cov))

            # Check if this is a duplicate.
            # The same conditions apply as in the old version.
            if 100 * idt > max_idt and \
                    100 * cov > max_aln_cov and \
                    abs(delta_l) < min_len_diff:
                is_duplicate = True
                log('    -> Duplicate!')
                break

        if is_duplicate == False:
            # This branch is not a duplicate. Add it to references,
            # so that the following branches can be compared to it afterwards.
            ref_seqs.append(compound_a_ctg[i])

            # Append the non-duplicates.
            new_header = ' '.join([a_ctg_id, v, w, len_, ovl, ne, str(delta_l), '%.2f' % (idt), '%.2f' % (cov)])
            ret.append((compound_a_ctg[i], new_header))

        log('')

    return ret

def run(max_idt, max_aln_cov, min_len_diff, ploidy):
    with open_fasta_reader("a_ctg_all.fa") as fp_in, \
         open("a_ctg.fa", "w") as fp_out:

        for compound_a_ctg in yield_single_compound(fp_in):
            filtered_a_ctg = filter_duplicate(compound_a_ctg, max_idt, max_aln_cov, min_len_diff, ploidy)

            for a_ctg, new_header in filtered_a_ctg:
                fp_out.write('>%s\n' % (new_header))
                fp_out.write(a_ctg.sequence)
                fp_out.write('\n')

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Removes duplicate a-tig, iff *all* conditions are violated. Assumes the working directory has the a_ctg_all.fa file, and produces a_ctg.fa',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--max-idt', type=int,
                        help="keep a-tig if the identity (in %) to the primary contig is <= max_idt", default=96)
    parser.add_argument('--max-aln-cov', type=int,
                        help="keep a-tig if the alignment coverage (in %) on the a-tig is <= max_aln_cov", default=97)
    parser.add_argument('--min-len-diff', type=int,
                        help="keep a-tig if the length different > min_len_diff", default=500)
    parser.add_argument('--ploidy', type=int,
                        help="For a diplid genome, 2 branches per SV are expected. This parameter limits the number of pairwise comparison. If <= 0, this threshold is not applied.", default=2)

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    run(args.max_idt, args.max_aln_cov, args.min_len_diff, args.ploidy)

if __name__ == "__main__":
    main(sys.argv)
