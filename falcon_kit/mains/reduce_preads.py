"""
Creates a reduced version of preads4falcon.fasta file by writing only the preads
which are incident with 'G' edges in the final assembly graph.
"""
from __future__ import absolute_import
from __future__ import print_function

import argparse
import sys
from falcon_kit.FastaReader import open_fasta_reader

def run(fp_out, pread_fasta_file, sg_edges_list_file):
    reads_in_layout = set()
    with open(sg_edges_list_file) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62 G"""
            v, w, rid, s, t, aln_score, idt, type_ = l
            if type_ != "G":
                continue
            r1 = v.split(":")[0]
            reads_in_layout.add(r1)
            r2 = w.split(":")[0]
            reads_in_layout.add(r2)

    with open_fasta_reader(pread_fasta_file) as f:
        for r in f:
            if r.name not in reads_in_layout:
                continue
            fp_out.write('>{}\n{}\n'.format(r.name, r.sequence.upper()))

def main(argv=sys.argv):
    description = 'Create a reduced set of preads, with only those used in the final layout.'
    parser = argparse.ArgumentParser(
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--pread_fasta_file', type=str,
            default='preads4falcon.fasta',
            help='Preads file, required to construct the contigs.')
    parser.add_argument('--sg_edges_list_file', type=str,
            default='sg_edges_list',
            help='File containing string graph edges, produced by ovlp_to_graph.py.')
    args = parser.parse_args(argv[1:])
    run(sys.stdout, **vars(args))

if __name__ == "__main__":
    main(sys.argv)
