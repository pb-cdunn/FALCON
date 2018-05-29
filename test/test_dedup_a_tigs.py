import falcon_kit.mains.dedup_a_tigs as mod
import helpers
import pytest
import os
import random

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_main_1(tmpdir, capsys):
    random.seed(1234567)

    dummy_seq_001_00 = ''.join([random.choice('ACTG') for i in xrange(45000)])
    dummy_seq_001_01 = dummy_seq_001_00[0:20000] + dummy_seq_001_00[25000:]

    dummy_seq_002_00 = dummy_seq_001_00
    dummy_seq_002_01 = dummy_seq_001_01
    dummy_seq_002_02 = dummy_seq_002_01

    fasta_lines = [
                    # This tests plain deduplication. The 001-01 is not a duplicate.
                    '>000000F-001-00 000000217:E 000000056:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq_001_00,
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq_001_01,
                    '>000000F-002-00 000000218:E 000000057:E 46974 1268418 33 0 1.00 1.00',
                    # This tests deduplication of 002-02 vs 002-01. The 002-02 should be removed.
                    dummy_seq_002_00,
                    '>000000F-002-01 000000124:E 000000079:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq_002_01,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq_002_02,
                  ]

    expected_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 -5000 1.00 0.50',
                    dummy_seq_001_01,
                    '>000000F-002-01 000000124:E 000000079:E 46974 1268418 33 -5000 1.00 0.50',
                    dummy_seq_002_01,
                  ]

    # Create a temporary a_ctg_all.fa file.
    test_out_a_ctg_all_file = tmpdir.join('a_ctg_all.fa')
    test_out_a_ctg_all_file.write('\n'.join(fasta_lines))

    argv = ['prog',
            '--a-ctg-all', str(test_out_a_ctg_all_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines) + '\n'

    assert(out == expected)

def test_main_2(tmpdir, capsys):
    """
    Test deduplication on very large sequences. This is supposed to test the
    heuristic of not-aligning very long sequences due to quadratic memory consumption
    of the current alignment procedure.
    Once the alignment approach is updated to support long sequences, this test
    should be updated to expect empty output (because the sequences are identical).
    """

    random.seed(1234567)

    # Generate a dummy sequence.
    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(300000)])
    fasta_lines = [
                    # Base sequence.
                    '>000000F-001-00 000000217:E 000000056:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    # Exactly the same sequence as a duplicate branch.
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq]

    # Create a temporary a_ctg_all.fa file.
    test_out_a_ctg_all_file = tmpdir.join('a_ctg_all.fa')
    test_out_a_ctg_all_file.write('\n'.join(fasta_lines))

    argv = ['prog',
            '--a-ctg-all', str(test_out_a_ctg_all_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    # The base sequence should never be output, and the "01" sequence should
    # ideally be deduplicated because it's identical to "00".
    # However, due to the quadratic memory in the current alignment implementation,
    # there is a max len threshold for alignment. If the sequences are longer
    # than this threshold, alignment and deduplication will be skipped, and the
    # associate contig output.
    # TODO: Once the alignment is updated to have linear memory, the `expected`
    # should be an empty string.
    expected = '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 -1.00 -1.00\n' + fasta_lines[3] + '\n'

    assert(out == expected)

def test_main_3(tmpdir, capsys):
    """
    Test deduplication on very short sequences.
    The legacy deduplication code had a 2000bp threshold on minimum sequence
    length for alignment. If any of the sequences was shorter than this,
    the alignment wouldn't be performed, and the sequence would not be deduplicated
    even if it was a duplicate.
    This tests for the same behaviour as the legacy code.
    The deduplication is designed as "conservative": keep the sequences unless
    there is strong evidence of duplication.
    TODO: Reconsider this logic.
    """

    random.seed(1234567)

    # Generate a dummy sequence.
    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(300)])
    fasta_lines = [
                    # Base sequence.
                    '>000000F-001-00 000000217:E 000000056:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    # Exactly the same sequence as a duplicate branch.
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq
                  ]

    expected_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 0.00 0.00',
                    dummy_seq,
                  ]

    # Create a temporary a_ctg_all.fa file.
    test_out_a_ctg_all_file = tmpdir.join('a_ctg_all.fa')
    test_out_a_ctg_all_file.write('\n'.join(fasta_lines))

    argv = ['prog',
            '--a-ctg-all', str(test_out_a_ctg_all_file),
            '--min-seq-len', '2000',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines) + '\n'

    assert(out == expected)

def test_main_4(tmpdir, capsys):
    """
    The DW alignment has a 100bp minimum distance threshold for (e1 - s1) and (e2 - s2).
    """

    random.seed(1234567)

    # Generate a dummy sequence.
    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(90)])
    fasta_lines = [
                    # Base sequence.
                    '>000000F-001-00 000000217:E 000000056:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    # Exactly the same sequence as a duplicate branch.
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq
                  ]

    expected_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 0.00 0.00',
                    dummy_seq,
                  ]

    # Create a temporary a_ctg_all.fa file.
    test_out_a_ctg_all_file = tmpdir.join('a_ctg_all.fa')
    test_out_a_ctg_all_file.write('\n'.join(fasta_lines))

    argv = ['prog',
            '--a-ctg-all', str(test_out_a_ctg_all_file),
            '--min-seq-len', '50',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines) + '\n'

    assert(out == expected)

def test_main_5(tmpdir, capsys):
    """
    Align two completely different sequences.
    """

    random.seed(1234567)

    # Generate a dummy sequence.
    dummy_seq_1 = ''.join([random.choice('AC') for i in xrange(3000)])
    dummy_seq_2 = ''.join([random.choice('GT') for i in xrange(3000)])
    fasta_lines = [
                    # Base sequence.
                    '>000000F-001-00 000000217:E 000000056:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq_1,
                    # Exactly the same sequence as a duplicate branch.
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq_2
                  ]

    expected_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 0.00 0.00',
                    dummy_seq_2,
                  ]

    # Create a temporary a_ctg_all.fa file.
    test_out_a_ctg_all_file = tmpdir.join('a_ctg_all.fa')
    test_out_a_ctg_all_file.write('\n'.join(fasta_lines))

    argv = ['prog',
            '--a-ctg-all', str(test_out_a_ctg_all_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines) + '\n'

    assert(out == expected)
