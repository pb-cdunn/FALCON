import falcon_kit.align_edlib as mod
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
    This works properly with Edlib. In the DW alignment module, this
    returns without alignment (since the sequences are too big).
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(300000)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 1.00, 1.00)  # delta_len, idt, cov

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



def test_get_global_aln_results_1():
    """
    Test aligning identical sequences.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(45000)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 1.00, 1.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_global_aln_results_2():
    """
    Test aligning non-identical sequences, where the query has a deletion
    compared to ref.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(45000)])
    query_seq = ref_seq[0:20000] + ref_seq[25000:]
    min_seq_len = 2000

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (-5000, 0.89, 1.0)  # delta_len, idt, cov

    assert(result == expected)

def test_get_global_aln_results_3():
    """
    Test aligning non-identical sequences, where the query has an insertion
    compared to ref.
    """

    query_seq = ''.join([random.choice('ACTG') for i in xrange(45000)])
    ref_seq = query_seq[0:20000] + query_seq[25000:]
    min_seq_len = 2000

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (5000, 0.89, 1.0)  # delta_len, idt, cov

    assert(result == expected)

def test_get_global_aln_results_4():
    """
    This tests alignment of very long sequences.
    Edlib is good with memory, and should align this easily.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(300000)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 1.00, 1.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_global_aln_results_5():
    """
    The legacy deduplication code has a threshold on minimum sequence
    length for alignment. If any of the sequences was shorter than this,
    the alignment wouldn't be performed, and the sequence would not be deduplicated
    even if it was a duplicate.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(300)])
    query_seq = ref_seq
    min_seq_len = 2000

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 0.00, 0.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_global_aln_results_6():
    """
    The DW alignment has a 100bp minimum distance threshold for (e1 - s1) and (e2 - s2).
    However, Edlib does not have any such constraints.
    """

    ref_seq = ''.join([random.choice('ACTG') for i in xrange(90)])
    query_seq = ref_seq
    min_seq_len = 50

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 1.00, 1.00)  # delta_len, idt, cov

    assert(result == expected)

def test_get_global_aln_results_7():
    """
    Align two completely different sequences.
    The align_edlib.py module will return coverage of 1.00 always, since
    the global alignment is applied.
    """

    ref_seq = ''.join([random.choice('AC') for i in xrange(3000)])
    query_seq = ''.join([random.choice('GT') for i in xrange(3000)])
    min_seq_len = 2000

    result = mod.get_global_aln_results(ref_seq, query_seq, min_seq_len)
    result = (result[0], round(result[1], 2), round(result[2], 2))
    expected = (0, 0.00, 1.00)  # delta_len, idt, cov

    assert(result == expected)

def test_count_cigar_ops_1():
    """
    Test empty input.
    """

    cigar = ''
    expected = (0, 0, 0, 0)           # num_m, num_i, num_d, total_len

    result = mod.count_cigar_ops(cigar)

    assert(result == expected)

def test_count_cigar_ops_2():
    """
    Test simple match ops, only one CIGAR op.
    """

    cigar = '10M'
    expected = (10, 0, 0, 10)           # num_m, num_i, num_d, total_len

    result = mod.count_cigar_ops(cigar)

    assert(result == expected)

def test_count_cigar_ops_3():
    """
    Test a more complex CIGAR string.
    """

    cigar = '10M1I123D4=2X'
    expected = (16, 1, 123, 140)           # num_m, num_i, num_d, total_len

    result = mod.count_cigar_ops(cigar)

    assert(result == expected)

def test_count_cigar_ops_4():
    """
    Test a degenerate case.
    """

    cigar = '10=X'
    expected = (0, 0, 0, 0)           # num_m, num_i, num_d, total_len

    with pytest.raises(Exception):
        result = mod.count_cigar_ops(cigar)

def test_count_cigar_ops_5():
    """
    Test a degenerate case.
    """

    cigar = '10I123'
    expected = (0, 0, 0, 0)           # num_m, num_i, num_d, total_len

    with pytest.raises(Exception):
        result = mod.count_cigar_ops(cigar)

def test_count_cigar_ops_6():
    """
    Test a degenerate case.
    """

    cigar = '12345'
    expected = (0, 0, 0, 0)           # num_m, num_i, num_d, total_len

    with pytest.raises(Exception):
        result = mod.count_cigar_ops(cigar)
