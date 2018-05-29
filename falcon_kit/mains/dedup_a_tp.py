from __future__ import absolute_import
from __future__ import print_function

from falcon_kit.FastaReader import open_fasta_reader
import argparse
import sys

def load_headers(fp_in):
    """
    Loads all a_ctg IDs from the a_ctg.fa, which is already deduplicated.
    """
    ret = set()
    for r in fp_in:
        a_ctg_id = r.name.split()[0]
        ret.add(a_ctg_id)
    return ret

def run(fp_out, a_ctg, a_ctg_all_tiling_path):
    with open_fasta_reader(a_ctg) as fp_in:
        a_ctg_ids = load_headers(fp_in)

    with open(a_ctg_all_tiling_path, 'r') as fp_in:
        for line in fp_in:
            line = line.strip()
            if len(line) == 0:  # pragma: no cover
                continue        # pragma: no cover
            sl = line.split()
            if sl[0] not in a_ctg_ids:
                continue
            fp_out.write('%s\n' % (line))

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Extracts all tiling paths from a_ctg_all_tiling_paths for which there is a header in a_ctg.fa (which was already deduplicated).',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--a-ctg', type=str,
                        help="Path to the a_ctg.fa file.", default='a_ctg.fa')
    parser.add_argument('--a-ctg-all-tiling-path', type=str,
                        help="Path to the a_ctg_all_tiling_path file.", default='a_ctg_all_tiling_path')

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    run(sys.stdout, **vars(args))

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
