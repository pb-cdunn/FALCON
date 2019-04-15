import falcon_kit.mains.fasta_filter as mod
import pytest
from falcon_kit import FastaReader
import functools
import helpers
import pytest
import os
import sys
import cStringIO
import json

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

#############################
### Define the test data. ###
#############################
fasta_tests = {}
expected_tests = {}
expected_general_median_tests = {}
expected_internal_median_tests = {}
expected_longest_tests = {}

fasta_tests[1] = ''
expected_tests[1] = ''

fasta_tests[2] = """\
>synthetic/1/0_500
GATTACA
>synthetic/2/0_500
GATTACA
>synthetic/3/0_500
GATTACA
>synthetic/4/0_500
GATTACA
>synthetic/5/0_500
GATTACA
"""
expected_tests[2] = fasta_tests[2]

fasta_tests[3] = """\
>synthetic/1/0_500
GATTACAGATTACA
>synthetic/1/0_500
GATTACAGATTACAGATTACAGATTACA
>synthetic/1/0_500
GATTACAGATTACAGATTACA
>synthetic/1/500_1000
GATTACAGATTACAGATTACAGATTACA
>synthetic/1/0_500
GATTACA
"""
expected_tests[3] = """\
>synthetic/1/0_500
GATTACAGATTACAGATTACA
"""

expected_longest_tests[3] = """\
>synthetic/1/0_500
GATTACAGATTACAGATTACAGATTACA
"""

fasta_tests[4] = """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/0_5
ACGTA
>synthetic/1/5_7
AC
>synthetic/1/7_10
ACG
>synthetic/1/10_14
ACGT
>synthetic/1/14_15
A
>synthetic/4/0_1
A
>synthetic/5/0_2
AC
>synthetic/5/2_5
ACG
>synthetic/5/5_9
ACGT
>synthetic/5/9_14
ACGTA
"""
expected_tests[4] = """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/7_10
ACG
>synthetic/4/0_1
A
>synthetic/5/5_9
ACGT
"""

fasta_tests[5] = """\
@synthetic/1/0_500
A
@synthetic/2/0_500
A
@synthetic/3/0_500
A
@synthetic/4/0_500
A
@synthetic/5/0_500'
A
"""
# That is invalid fasta.
expected_tests[5] = None    # The module under test should throw on the above input

fasta_tests[6] = """\
>synthetic/5/2_5
ACG
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/0_5
ACGTA
>synthetic/1/10_14
ACGT
>synthetic/1/14_15
A
>synthetic/4/0_1
A
>synthetic/5/0_2
AC
>synthetic/5/5_9
ACGT
>synthetic/5/9_14
ACGTA
>synthetic/1/5_7
AC
>synthetic/1/7_10
ACG
"""
expected_tests[6] = """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/4/0_1
A
>synthetic/5/5_9
ACGT
>synthetic/1/7_10
ACG
"""

# Test ccs input
fasta_tests[8] = """\
>synthetic/1/ccs
GATTACA
"""
expected_tests[8] = """\
>synthetic/1/0_7
GATTACA
"""

# Specific for the internal-median command.
fasta_tests[7] = """\
>synthetic/1/0_3
AAA
>synthetic/1/3_7
ACAC
>synthetic/1/7_10
AAA
>synthetic/2/0_3
TTT
>synthetic/2/3_8
TTTGG
>synthetic/3/0_5
TTTGG
>synthetic/3/5_8
TTT
>synthetic/4/0_5
CCCCC
"""
expected_general_median_tests[7] = """\
>synthetic/1/7_10
AAA
>synthetic/2/3_8
TTTGG
>synthetic/3/0_5
TTTGG
>synthetic/4/0_5
CCCCC
"""
expected_internal_median_tests[7] = """\
>synthetic/1/3_7
ACAC
>synthetic/2/3_8
TTTGG
>synthetic/3/0_5
TTTGG
>synthetic/4/0_5
CCCCC
"""

whitelist_tests_general_median = [
    {   'fasta': fasta_tests[4],
        'whitelist': set(),
        'expected': expected_tests[4]
    },
    {   'fasta': fasta_tests[4],
        'whitelist': None,
        'expected': expected_tests[4]
    },
    {   'fasta': fasta_tests[4],
        'whitelist': set(['synthetic/2', 'synthetic/4', 'synthetic/3']),
        'expected': """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/4/0_1
A
"""
    },
    {   'fasta': fasta_tests[4],
        'whitelist': set(['synthetic/2', 'synthetic/3', 'synthetic/1', 'synthetic/4', 'synthetic/5']),
        'expected': expected_tests[4]
    },
]

whitelist_tests_longest = [
    {   'fasta': fasta_tests[4],
        'whitelist': set(),
        'expected': """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/0_5
ACGTA
>synthetic/4/0_1
A
>synthetic/5/9_14
ACGTA
"""
    },
    {   'fasta': fasta_tests[4],
        'whitelist': None,
        'expected': """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/0_5
ACGTA
>synthetic/4/0_1
A
>synthetic/5/9_14
ACGTA
"""
    },
    {   'fasta': fasta_tests[4],
        'whitelist': set(['synthetic/2', 'synthetic/4', 'synthetic/3']),
        'expected': """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/4/0_1
A
"""
    },
    {   'fasta': fasta_tests[4],
        'whitelist': set(['synthetic/2', 'synthetic/3', 'synthetic/1', 'synthetic/4', 'synthetic/5']),
        'expected': """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/0_5
ACGTA
>synthetic/4/0_1
A
>synthetic/5/9_14
ACGTA
"""
    },
]

whitelist_tests_pass = [
    {   'fasta': fasta_tests[4],
        'whitelist': set(),
        'expected': fasta_tests[4]
    },
    {   'fasta': fasta_tests[4],
        'whitelist': None,
        'expected': fasta_tests[4]
    },
    {   'fasta': fasta_tests[4],
        'whitelist': set(['synthetic/2', 'synthetic/1', 'synthetic/3']),
        'expected': """\
>synthetic/2/0_1
A
>synthetic/3/0_1
A
>synthetic/1/0_5
ACGTA
>synthetic/1/5_7
AC
>synthetic/1/7_10
ACG
>synthetic/1/10_14
ACGT
>synthetic/1/14_15
A
"""
    },
    {   'fasta': fasta_tests[4],
        'whitelist': set(['synthetic/2', 'synthetic/3', 'synthetic/1', 'synthetic/4', 'synthetic/5']),
        'expected': fasta_tests[4]
    },
]

#########################
### Helper functions. ###
#########################
def check_run(func, fasta, expected, whitelist_set=None):
    if expected is None:
        return check_run_raises(func, fasta, whitelist_set)
    fp_in = cStringIO.StringIO(fasta)
    fp_out = cStringIO.StringIO()
    func(fp_in, fp_out, whitelist_set)
    if expected != fp_out.getvalue():
        print(fp_out.getvalue())
    assert expected == fp_out.getvalue()

def check_run_raises(func, fasta, whitelist_set=None):
    fp_in = cStringIO.StringIO(fasta)
    fp_out = cStringIO.StringIO()
    with pytest.raises(Exception) as e:
        func(fp_in, fp_out, whitelist_set)
    assert not isinstance(e.value, AssertionError) # Test assertions elsewhere.

def check_main_from_file(cmd, fasta, expected, in_fa_file, capsys):
    in_fa_file.write(fasta)
    argv = ['prog', cmd, str(in_fa_file)]
    mod.main(argv)
    out, err = capsys.readouterr()
    if expected != out:
        print(out)
    assert(out == expected)

def check_main_from_file_with_whitelist(cmd, fasta, expected, in_whitelist_file, in_fa_file, capsys):
    in_fa_file.write(fasta)
    argv = ['prog', '--zmw-whitelist-fn', in_whitelist_file, cmd, str(in_fa_file)]
    mod.main(argv)
    out, err = capsys.readouterr()
    if expected != out:
        print(out)
    assert(out == expected)

def check_main_from_stdin(cmd, fasta, expected, capsys):
    mod.sys.stdin = cStringIO.StringIO(fasta)
    argv = ['prog', cmd, '-']
    mod.main(argv)
    out, err = capsys.readouterr()
    if expected != out:
        print(out)
    assert(out == expected)

def check_main_from_stdin_with_whitelist(cmd, fasta, expected, in_whitelist_file, capsys):
    mod.sys.stdin = cStringIO.StringIO(fasta)
    argv = ['prog', '--zmw-whitelist-fn',  in_whitelist_file, cmd, '-']
    mod.main(argv)
    out, err = capsys.readouterr()
    if expected != out:
        print(out)
    assert(out == expected)

#############################
## Test internal functions ##
#############################
def test_yield_record_and_tokenized_headers():
    # We could delete this test, since we have whitelist tests also.
    # Or we could expand this.
    fp_in = cStringIO.StringIO(fasta_tests[1])
    records = FastaReader.yield_fasta_record(fp_in, log=lambda x: None)
    w = set()
    assert [] == list(mod.yield_record_and_tokenized_headers(w, records))

def test_load_zmw_whitelist(tmpdir, caplog):
    fn = ''
    expected = set()
    assert mod.load_zmw_whitelist(fn) == expected
    assert not caplog.records

    path = tmpdir.join('whitelist_empty.json')
    path.write('') # empty
    fn = str(path)
    expected = set()
    assert mod.load_zmw_whitelist(fn) == set()
    assert 'Assuming empty whitelist' in caplog.text
    caplog.clear()

    path = tmpdir.join('whitelist_empty_list.json')
    path.write('[]') # empty
    fn = str(path)
    expected = set()
    assert mod.load_zmw_whitelist(fn) == expected
    assert not caplog.records

    redundant = ['aa/1', 'bb/2', 'aa/1']
    path = tmpdir.join('whitelist_redundant_list.json')
    path.write(json.dumps(redundant))
    fn = str(path)
    expected = set(['aa/1', 'bb/2'])
    assert mod.load_zmw_whitelist(fn) == expected
    assert not caplog.records

#############################
### Test the PASS filter. ###
#############################
@pytest.mark.parametrize('case', [
    1, 2, 3, 4,
])
def test_run_pass_filter(case):
    provided = expected = fasta_tests[case]
    check_run(mod.run_pass_filter, provided, expected)

def test_run_pass_filter_5():
    check_run_raises(mod.run_pass_filter, fasta_tests[5])

def test_run_pass_filter_8_ccs():
    provided = fasta_tests[8]
    expected = expected_tests[8]
    check_run(mod.run_pass_filter, provided, expected)

########################################
### Test the streamed median filter. ###
########################################
@pytest.mark.parametrize('case', [
    1, 2, 3, 4, 5,
])
def test_run_streamed_median_filter(case):
    provided = fasta_tests[case]
    expected = expected_tests[case]
    check_run(mod.run_streamed_median_filter, provided, expected)

def test_run_streamed_median_7g():
    # test general
    curried = functools.partial(mod.run_streamed_median_filter, zmw_filter_func=mod.median_zmw_subread)
    check_run(curried, fasta_tests[7], expected_general_median_tests[7])

def test_run_streamed_median_7i():
    # test internal
    curried = functools.partial(mod.run_streamed_median_filter, zmw_filter_func=mod.internal_median_zmw_subread)
    check_run(curried, fasta_tests[7], expected_internal_median_tests[7])

def test_run_streamed_longest_3():
    check_run(mod.run_streamed_longest_filter, fasta_tests[3], expected_longest_tests[3])

#######################################
### Test the general median filter. ###
#######################################
@pytest.mark.parametrize('case', [
    1, 2, 3, 4, 5, 6
])
def test_run_median_filter(case):
    provided = fasta_tests[case]
    expected = expected_tests[case]
    check_run(mod.run_median_filter, provided, expected)

def test_run_median_filter_piped_input():
    """
    The run_median_filter asserts if the input file doesn't exist.
    This is to prevent specifying streams we can't rewind.
    """
    with pytest.raises(AssertionError) as e:
        mod.run_median_filter(sys.stdin, '', '')
    assert 'Cannot rewind' in str(e.value)

def test_run_median_filter_7():
    check_run(mod.run_median_filter, fasta_tests[7], expected_general_median_tests[7])

########################################
### Test the internal median filter. ###
########################################
@pytest.mark.parametrize('case', [
    1, 2, 3, 4, 5, 6
])
def test_run_internal_median_filter(case):
    provided = fasta_tests[case]
    expected = expected_tests[case]
    check_run(mod.run_internal_median_filter, provided, expected)

def test_run_internal_median_filter_piped_input():
    """
    The run_median_filter asserts if the input file doesn't exist.
    This is to prevent specifying streams we can't rewind.
    """
    check_run_raises(mod.run_internal_median_filter, fasta_tests[5])
    with pytest.raises(AssertionError) as e:
        mod.run_internal_median_filter(sys.stdin, '', '')
    assert 'Cannot rewind' in str(e.value)

def test_run_internal_median_filter_7():
    check_run(mod.run_internal_median_filter, fasta_tests[7], expected_internal_median_tests[7])

###############################
### Test whitelisting ZMWs. ###
###############################
@pytest.mark.parametrize('test_data', whitelist_tests_general_median)
def test_whitelist_run_internal_median_filter(test_data):
    test_fasta = test_data['fasta']
    test_whitelist = test_data['whitelist']
    test_expected = test_data['expected']
    check_run(mod.run_internal_median_filter, test_fasta, test_expected, test_whitelist)

@pytest.mark.parametrize('test_data', whitelist_tests_general_median)
def test_whitelist_run_median_filter(test_data):
    test_fasta = test_data['fasta']
    test_whitelist = test_data['whitelist']
    test_expected = test_data['expected']
    check_run(mod.run_median_filter, test_fasta, test_expected, test_whitelist)

@pytest.mark.parametrize('test_data', whitelist_tests_general_median)
def test_whitelist_run_streamed_median_filter(test_data):
    test_fasta = test_data['fasta']
    test_whitelist = test_data['whitelist']
    test_expected = test_data['expected']
    check_run(mod.run_streamed_median_filter, test_fasta, test_expected, test_whitelist)

@pytest.mark.parametrize('test_data', whitelist_tests_longest)
def test_whitelist_run_streamed_longest_filter(test_data):
    test_fasta = test_data['fasta']
    test_whitelist = test_data['whitelist']
    test_expected = test_data['expected']
    check_run(mod.run_streamed_longest_filter, test_fasta, test_expected, test_whitelist)

@pytest.mark.parametrize('test_data', whitelist_tests_pass)
def test_whitelist_run_pass_filter(test_data):
    test_fasta = test_data['fasta']
    test_whitelist = test_data['whitelist']
    test_expected = test_data['expected']
    check_run(mod.run_pass_filter, test_fasta, test_expected, test_whitelist)

###############################
### Test the main commands. ###
###############################
def test_main_cmd_pass_1(capsys):
    check_main_from_stdin('pass', fasta_tests[4], fasta_tests[4], capsys)

def test_main_cmd_pass_2(tmpdir, capsys):
    in_fa_file = tmpdir.join('in.fa')
    check_main_from_file('pass', fasta_tests[4], fasta_tests[4], in_fa_file, capsys)

def test_main_cmd_streamed_median_1(capsys):
    check_main_from_stdin('streamed-median', fasta_tests[4], expected_tests[4], capsys)

def test_main_cmd_streamed_median_2(tmpdir, capsys):
    in_fa_file = tmpdir.join('in.fa')
    check_main_from_file('streamed-median', fasta_tests[4], expected_tests[4], in_fa_file, capsys)

def test_main_cmd_median_1(tmpdir, capsys):
    in_fa_file = tmpdir.join('in.fa')
    check_main_from_file('median', fasta_tests[4], expected_tests[4], in_fa_file, capsys)

def test_main_cmd_internal_median_1(tmpdir, capsys):
    in_fa_file = tmpdir.join('in.fa')
    check_main_from_file('internal-median', fasta_tests[4], expected_tests[4], in_fa_file, capsys)

def test_main_cmd_streamed_internal_median_1(tmpdir, capsys):
    check_main_from_stdin('streamed-internal-median', fasta_tests[4], expected_tests[4], capsys)

def test_main_cmd_streamed_internal_median_2(tmpdir, capsys):
    in_fa_file = tmpdir.join('in.fa')
    check_main_from_file('streamed-internal-median', fasta_tests[4], expected_tests[4], in_fa_file, capsys)

#####################################################
### Test ZMW whitelisting from the main commands. ###
#####################################################
def test_whitelist_main_cmd_streamed_longest_from_stdin(tmpdir, capsys):
    for test_id, test_data in enumerate(whitelist_tests_longest):
        test_fasta = test_data['fasta']
        test_whitelist = test_data['whitelist']
        test_expected = test_data['expected']

        test_whitelist_fn = tmpdir.join('zmw.whitelist.{}.json'.format(test_id))
        if test_whitelist == None:
            # In case None, just omit the path.
            test_whitelist_fn = ''
        else:
            # In case an empty list, write an empty file.
            test_whitelist_fn.write(json.dumps(list(test_whitelist)))

        check_main_from_stdin_with_whitelist('streamed-longest', test_fasta, test_expected, str(test_whitelist_fn), capsys)

#######################
### Unit-test functions
def test_tokenize_header():
    with pytest.raises(ValueError) as e:
        mod.tokenize_header('bad_header')
    assert 'bad_header' in str(e.value)
