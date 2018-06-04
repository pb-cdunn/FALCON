from __future__ import absolute_import
from __future__ import print_function

from falcon_kit.FastaReader import open_fasta_reader
import argparse
import sys
# import falcon_kit.align_dw as align
import falcon_kit.align_edlib as align

class TooLongError(Exception): pass

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')

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

def filter_duplicate(compound_a_ctg, max_idt, max_aln_cov, min_len_diff, min_seq_len, ploidy):
    """
    Takes a list of a_ctg sequences in a compound unitig (bubble) which need
    to be deduplicated according to the parameters.
    The zeroth sequence in the list is the "base" sequence. This sequence is
    already part of the primary path by definition, and will not be output.
    """

    ret = []

    # Sanity check.
    if len(compound_a_ctg) == 0: return ret # pragma: no cover

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
            delta_l, idt, cov = align.get_aln_results(ref_seq, a_ctg_seq, min_seq_len)

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

def run(fp_out, fp_in, max_idt, max_aln_cov, min_len_diff, min_seq_len, ploidy):
    for compound_a_ctg in yield_single_compound(fp_in):
        filtered_a_ctg = filter_duplicate(compound_a_ctg, max_idt, max_aln_cov, min_len_diff, min_seq_len, ploidy)

        for a_ctg, new_header in filtered_a_ctg:
            fp_out.write('>%s\n' % (new_header))
            fp_out.write(a_ctg.sequence)
            fp_out.write('\n')

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Removes duplicate a-tig, iff *all* conditions are violated. Assumes the working directory has the a_ctg_all.fa file, and produces a_ctg.fa',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--max-idt', type=int,
                        help="Keep a-tig if the identity (in %%) to the primary contig is <= max_idt", default=96)
    parser.add_argument('--max-aln-cov', type=int,
                        help="Keep a-tig if the alignment coverage (in %%) on the a-tig is <= max_aln_cov", default=97)
    parser.add_argument('--min-len-diff', type=int,
                        help="Keep a-tig if the length different > min_len_diff", default=500)
    parser.add_argument('--min-seq-len', type=int,
                        help="Branches with length less than this threshold will always be deduplicated.", default=2000)
    parser.add_argument('--ploidy', type=int,
                        help="For a diplid genome, 2 branches per SV are expected. This parameter limits the number of pairwise comparison. If <= 0, this threshold is not applied.", default=2)
    parser.add_argument('--a-ctg-all', type=str,
                        help="Input set of all associate contigs for deduplication.", default="a_ctg_all.fa")

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    with open_fasta_reader(args.a_ctg_all) as fp_in:
        run(sys.stdout, fp_in, args.max_idt, args.max_aln_cov, args.min_len_diff, args.min_seq_len, args.ploidy)

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
