"""
Creates a reduced version of preads4falcon.fasta file by writing only the preads
which are incident with 'G' edges in the final assembly graph.
"""
from __future__ import absolute_import
from __future__ import print_function

import argparse
import logging
import sys
from ..FastaReader import open_fasta_reader
from ..io import open_progress

default_sg_edges_list_fns = ['./sg_edges_list']

def run(fp_out, preads_fasta_fn, sg_edges_list_fns):
    # Workaround the Argparse issue. It does not override
    # the default argument value when the parameter is
    # used in the append mode, but instead adds to the default
    # list. https://bugs.python.org/issue16399
    # Instead, we will not specify the default value, and
    # check if the list is emptu here here, so that the user
    # can specify exactly the paths to the file(s).
    if not sg_edges_list_fns:
        sg_edges_list_fns = default_sg_edges_list_fns

    reads_in_layout = set()

    for fn in sg_edges_list_fns:
        with open_progress(fn) as fp_in:
            for l in fp_in:
                l = l.strip().split()
                """001039799:E 000333411:E 000333411 17524 20167 17524 99.62 G"""
                v, w, rid, s, t, aln_score, idt, type_ = l
                if type_ != "G":
                    continue
                r1 = v.split(":")[0]
                reads_in_layout.add(r1)
                r2 = w.split(":")[0]
                reads_in_layout.add(r2)

    with open_fasta_reader(preads_fasta_fn) as f:
        for r in f:
            if r.name not in reads_in_layout:
                continue
            fp_out.write('>{}\n{}\n'.format(r.name, r.sequence.upper()))

def main(argv=sys.argv):
    description = 'Create a reduced set of preads, with only those used in the final layout. Write to stdout.'
    parser = argparse.ArgumentParser(
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--preads-fasta-fn', type=str,
            default='preads4falcon.fasta',
            help='Preads file, required to construct the contigs.')
    parser.add_argument('--sg-edges-list-fns', action='append',
            help='One or more files containing string graph edges, produced by ovlp_to_graph.py.')
    args = parser.parse_args(argv[1:])
    run(sys.stdout, **vars(args))

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main(sys.argv)
