import falcon_kit.mains.dedup_a_tp as mod
import helpers
import pytest
import os
import random
from falcon_kit.FastaReader import open_fasta_reader

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_main_1(tmpdir, capsys):
    """
    Empty input files should produce empty output.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                   ]

    tp_lines = [
                ]

    expected_lines = [
                    ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    test_out_a_ctg_all_tp_file = tmpdir.join('a_ctg_all_tiling_path')
    test_out_a_ctg_all_tp_file.write('\n'.join(tp_lines) + '\n')
    # test_out_a_ctg_all_tp_file.write('\n')

    argv = ['prog',
            '--a-ctg', str(test_out_a_ctg_file),
            '--a-ctg-all-tiling-path', str(test_out_a_ctg_all_tp_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines)

    assert(out == expected)

def test_main_2(tmpdir, capsys):
    """
    Regular test case.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                  ]

    tp_lines = [
                    '000000F-001-00 000000217:E 000000189:E 000000189 33799 36467 33798 99.97',
                    '000000F-001-00 000000189:E 000000223:E 000000223 35637 39957 35637 99.96',
                    '000000F-001-00 000000223:E 000000273:E 000000273 38476 39968 38476 99.96',

                    '000000F-001-01 000000217:E 000000110:E 000000110 33829 37916 33802 99.86',
                    '000000F-001-01 000000110:E 000000037:E 000000037 35944 39971 35947 99.97',
                    '000000F-001-01 000000037:E 000000131:E 000000131 39313 39968 39313 99.97',

                    '000000F-002-00 000000035:E 000000078:E 000000078 38260 39977 38260 99.98',
                    '000000F-002-00 000000078:E 000000020:E 000000020 38490 39969 38490 99.97',
                    '000000F-002-00 000000020:E 000000142:E 000000142 39547 39983 39547 99.98',

                    '000000F-002-01 000000035:E 000000228:E 000000228 38670 39961 38670 99.89',
                    '000000F-002-01 000000228:E 000000201:E 000000201 38578 39971 38571 99.95',
                    '000000F-002-01 000000201:E 000000243:E 000000243 39822 39967 39822 99.96',

                    '000000F-002-02 000000176:E 000000250:E 000000250 39602 39966 39602 99.97',
                    '000000F-002-02 000000250:E 000000211:E 000000211 38613 39982 38610 99.98',
                    '000000F-002-02 000000211:E 000000177:E 000000177 39666 39984 39668 99.98',
                ]

    expected_lines = [
                    '000000F-001-01 000000217:E 000000110:E 000000110 33829 37916 33802 99.86',
                    '000000F-001-01 000000110:E 000000037:E 000000037 35944 39971 35947 99.97',
                    '000000F-001-01 000000037:E 000000131:E 000000131 39313 39968 39313 99.97',

                    '000000F-002-02 000000176:E 000000250:E 000000250 39602 39966 39602 99.97',
                    '000000F-002-02 000000250:E 000000211:E 000000211 38613 39982 38610 99.98',
                    '000000F-002-02 000000211:E 000000177:E 000000177 39666 39984 39668 99.98',
                ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    test_out_a_ctg_all_tp_file = tmpdir.join('a_ctg_all_tiling_path')
    test_out_a_ctg_all_tp_file.write('\n'.join(tp_lines) + '\n')
    # test_out_a_ctg_all_tp_file.write('\n')

    argv = ['prog',
            '--a-ctg', str(test_out_a_ctg_file),
            '--a-ctg-all-tiling-path', str(test_out_a_ctg_all_tp_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines) + '\n'

    assert(out == expected)


def test_main_3(tmpdir, capsys):
    """
    Empty tiling path file.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                  ]

    tp_lines = [
                ]

    expected_lines = [
                    ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    test_out_a_ctg_all_tp_file = tmpdir.join('a_ctg_all_tiling_path')
    test_out_a_ctg_all_tp_file.write('\n'.join(tp_lines) + '\n')
    # test_out_a_ctg_all_tp_file.write('\n')

    argv = ['prog',
            '--a-ctg', str(test_out_a_ctg_file),
            '--a-ctg-all-tiling-path', str(test_out_a_ctg_all_tp_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines)

    assert(out == expected)

def test_main_4(tmpdir, capsys):
    """
    Tiling path does not contain all contigs present in a_ctg.fa. Only the ones
    which match should be output, still.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                  ]

    tp_lines = [
                    '000000F-001-01 000000217:E 000000110:E 000000110 33829 37916 33802 99.86',
                    '000000F-001-01 000000110:E 000000037:E 000000037 35944 39971 35947 99.97',
                    '000000F-001-01 000000037:E 000000131:E 000000131 39313 39968 39313 99.97',
                ]

    expected_lines = [
                    '000000F-001-01 000000217:E 000000110:E 000000110 33829 37916 33802 99.86',
                    '000000F-001-01 000000110:E 000000037:E 000000037 35944 39971 35947 99.97',
                    '000000F-001-01 000000037:E 000000131:E 000000131 39313 39968 39313 99.97',
                ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    test_out_a_ctg_all_tp_file = tmpdir.join('a_ctg_all_tiling_path')
    test_out_a_ctg_all_tp_file.write('\n'.join(tp_lines) + '\n')
    # test_out_a_ctg_all_tp_file.write('\n')

    argv = ['prog',
            '--a-ctg', str(test_out_a_ctg_file),
            '--a-ctg-all-tiling-path', str(test_out_a_ctg_all_tp_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = '\n'.join(expected_lines) + '\n'

    assert(out == expected)

def test_load_headers_1(tmpdir):
    """
    Regular case.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                  ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    # fp_in = StringIO('\n'.join(fasta_lines))
    with open_fasta_reader(str(test_out_a_ctg_file)) as fp_in:
        out = mod.load_headers(fp_in)

    expected = set(['000000F-001-01', '000000F-002-02'])

    assert(out == expected)

def test_load_headers_2(tmpdir):
    """
    There is a duplicate header in the input. This should be reduced to
    only the unique occurrences.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                    '>000000F-001-01 000000123:E 000000078:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                    '>000000F-002-02 000000125:E 000000080:E 46974 1268418 33 0 1.00 1.00',
                    dummy_seq,
                  ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    # fp_in = StringIO('\n'.join(fasta_lines))
    with open_fasta_reader(str(test_out_a_ctg_file)) as fp_in:
        out = mod.load_headers(fp_in)

    expected = set(['000000F-001-01', '000000F-002-02'])

    assert(out == expected)

def test_load_headers_3(tmpdir):
    """
    Empty input.
    """
    random.seed(1234567)

    dummy_seq = ''.join([random.choice('ACTG') for i in xrange(500)])

    fasta_lines = [
                  ]

    # Create a temporary input files.
    test_out_a_ctg_file = tmpdir.join('a_ctg.fa')
    test_out_a_ctg_file.write('\n'.join(fasta_lines))

    # fp_in = StringIO('\n'.join(fasta_lines))
    with open_fasta_reader(str(test_out_a_ctg_file)) as fp_in:
        out = mod.load_headers(fp_in)

    expected = set()

    assert(out == expected)
