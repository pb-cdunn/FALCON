import falcon_kit.align_dw as mod
import helpers
import pytest
import os
import random

random.seed(1234567)

def test_get_aln_results_1():
    """
    Test aligning identical sequences.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(45000)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 1.00, 1.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_aln_results_2():
    """
    Test aligning non-identical sequences.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(45000)])
    query_seq = ref_seq[0:20000] + ref_seq[25000:]
    min_seq_len = 2000

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (-5000, 1.00, 0.50)  # delta_len, idt, cov

    assert(result == expected)

def test_get_aln_results_3():
    """
    This tests alignment of very long sequences.
    Ideally, this would return (0, 1.00, 1.00), but the DW alignment
    implementation has a quadratic memory complexity, so there is a
    hardcoded threshold than either ref or query cannot be longer than
    250000 bp. In this case, (0, -1.0, -1.0) is returned.
    This test case would fail with a normal linear-memory aligner,
    which is actually what would be expected.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(300000)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, -1.00, -1.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_aln_results_4():
    """
    The legacy deduplication code has a threshold on minimum sequence
    length for alignment. If any of the sequences was shorter than this,
    the alignment wouldn't be performed, and the sequence would not be deduplicated
    even if it was a duplicate.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(300)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 0.00, 0.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_aln_results_5():
    """
    The DW alignment has a 100bp minimum distance threshold for (e1 - s1) and (e2 - s2).
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(90)])
    query_seq = ref_seq
    min_seq_len = 50

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 0.00, 0.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_aln_results_6():
    """
    Align two completely different sequences.
    """

    ref_seq = ''.join([random.choice('AC') for i in xrange(3000)])
    query_seq = ''.join([random.choice('GT') for i in xrange(3000)])
    min_seq_len = 2000

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 0.00, 0.00)  # delta_len, idt, cov

    assert(result == expected)
