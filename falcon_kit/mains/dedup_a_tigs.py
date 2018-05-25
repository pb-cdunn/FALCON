from __future__ import absolute_import
from __future__ import print_function

from falcon_kit.FastaReader import open_fasta_reader
import argparse
import sys
from falcon_kit import kup, falcon, DWA

def get_aln_results(ref_seq, query_seq):
    # Align the a_ctg against the base.
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
            log('WARNING: Seqs were too long for get_aln_data(), so we set idt/cov low enough to prevent filtering by dedup_a_tigs, at atig_path[:-1] == {}'.format(atig_path[:-1]))
            idt = -1.0
            cov = -1.0
    return delta_len, idt, cov

def run(max_idt, max_aln_cov, min_len_diff):
    with open_fasta_reader("a_ctg_all.fa") as reads:
        with open("a_ctg.fa", "w") as f:
            for r in reads:
                tig_id, v, w, len_, ovl, ne, delta_l, idt, cov = r.name.split()
                if 100 * float(idt) > max_idt and 100 * float(cov) > max_aln_cov and\
                   abs(int(delta_l)) < min_len_diff:
                    continue
                print(">" + r.name, file=f)
                print(r.sequence, file=f)

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Removes duplicate a-tig, iff *all* conditions are violated. Assumes the working directory has the a_ctg_all.fa file, and produces a_ctg.fa',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--max-idt', type=int,
                        help="keep a-tig if the identity (in %) to the primary contig is <= max_idt", default=96)
    parser.add_argument('--max-aln-cov', type=int,
                        help="keep a-tig if the alignment coverage (in %) on the a-tig is <= max_aln_cov", default=97)
    parser.add_argument('--min-len-diff', type=int,
                        help="keep a-tig if the length different > min_len_diff", default=500)
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))

if __name__ == "__main__":
    main(sys.argv)
