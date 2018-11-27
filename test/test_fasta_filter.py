import falcon_kit.mains.fasta_filter as mod
import functools
import helpers
import pytest
import os
from falcon_kit.FastaReader import open_fasta_reader
import cStringIO

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
GATTACAGATTACAGATTACAGATTACA
>synthetic/1/0_500
GATTACAGATTACA
>synthetic/1/0_500
GATTACAGATTACAGATTACA
>synthetic/1/0_500
GATTACAGATTACAGATTACAGATTACA
>synthetic/1/0_500
GATTACA
"""
expected_tests[3] = """\
>synthetic/1/0_500
GATTACAGATTACAGATTACA
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

#########################
### Helper functions. ###
#########################
def check_run(func, fasta, expected, fasta_path):
    fp_in = cStringIO.StringIO(fasta)
    fp_out = cStringIO.StringIO()
    func(fp_in, fp_out, fasta_path)
    if expected != fp_out.getvalue():
        print(fp_out.getvalue())
    assert expected == fp_out.getvalue()

def check_run_raises(func, fasta, fasta_path):
    fp_in = cStringIO.StringIO(fasta)
    with pytest.raises(Exception):
        func(fp_in, fp_out, fasta_path)

def check_main_from_file(cmd, fasta, expected, in_fa_file, capsys):
    in_fa_file.write(fasta)
    argv = ['prog', cmd, str(in_fa_file)]
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

#############################
### Test the PASS filter. ###
#############################
def test_run_pass_filter_1():
    check_run(mod.run_pass_filter, fasta_tests[1], fasta_tests[1], '-')

def test_run_pass_filter_2():
    check_run(mod.run_pass_filter, fasta_tests[2], fasta_tests[2], '-')

def test_run_pass_filter_3():
    check_run(mod.run_pass_filter, fasta_tests[3], fasta_tests[3], '-')

def test_run_pass_filter_4():
    check_run(mod.run_pass_filter, fasta_tests[4], fasta_tests[4], '-')

def test_run_pass_filter_5():
    check_run_raises(mod.run_pass_filter, fasta_tests[5], '-')

########################################
### Test the streamed median filter. ###
########################################
def test_run_streamed_median_1():
    check_run(mod.run_streamed_median_filter, fasta_tests[1], expected_tests[1], '-')

def test_run_streamed_median_2():
    check_run(mod.run_streamed_median_filter, fasta_tests[2], expected_tests[2], '-')

def test_run_streamed_median_3():
    check_run(mod.run_streamed_median_filter, fasta_tests[3], expected_tests[3], '-')

def test_run_streamed_median_4():
    check_run(mod.run_streamed_median_filter, fasta_tests[4], expected_tests[4], '-')

def test_run_streamed_median_5():
    check_run_raises(mod.run_streamed_median_filter, fasta_tests[5], '-')

def test_run_streamed_median_7g():
    # test general
    curried = functools.partial(mod.run_streamed_median_filter, zmw_filter_func=mod.median_zmw_subread)
    check_run(curried, fasta_tests[7], expected_general_median_tests[7], '-')

def test_run_streamed_median_7i():
    # test internal
    curried = functools.partial(mod.run_streamed_median_filter, zmw_filter_func=mod.internal_median_zmw_subread)
    check_run(curried, fasta_tests[7], expected_internal_median_tests[7], '-')

#######################################
### Test the general median filter. ###
#######################################
def test_run_median_filter_1(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[1]
    in_fa_file.write(fasta)
    check_run(mod.run_median_filter, fasta, expected_tests[1], str(in_fa_file))

def test_run_median_filter_2(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[2]
    in_fa_file.write(fasta)
    check_run(mod.run_median_filter, fasta, expected_tests[2], str(in_fa_file))

def test_run_median_filter_3(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[3]
    in_fa_file.write(fasta)
    check_run(mod.run_median_filter, fasta, expected_tests[3], str(in_fa_file))

def test_run_median_filter_4(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[4]
    in_fa_file.write(fasta)
    check_run(mod.run_median_filter, fasta, expected_tests[4], str(in_fa_file))

def test_run_median_filter_5(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[5]
    in_fa_file.write(fasta)
    check_run_raises(mod.run_median_filter, fasta_tests[5], str(in_fa_file))

def test_run_median_filter_6(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[6]
    in_fa_file.write(fasta)
    check_run(mod.run_median_filter, fasta_tests[6], expected_tests[6], str(in_fa_file))

def test_run_median_filter_5x(tmpdir):
    """
    The run_median_filter asserts if the input file doesn't exist.
    This is to prevent specifying streams we can't rewind.
    """
    check_run_raises(mod.run_median_filter, fasta_tests[5], '-')

def test_run_median_filter_7(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[6]
    in_fa_file.write(fasta)
    check_run(mod.run_median_filter, fasta_tests[7], expected_general_median_tests[7], str(in_fa_file))

########################################
### Test the internal median filter. ###
########################################
def test_run_internal_median_filter_1(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[1]
    in_fa_file.write(fasta)
    check_run(mod.run_internal_median_filter, fasta, expected_tests[1], str(in_fa_file))

def test_run_minternal_edian_filter_2(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[2]
    in_fa_file.write(fasta)
    check_run(mod.run_internal_median_filter, fasta, expected_tests[2], str(in_fa_file))

def test_run_internal_median_filter_3(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[3]
    in_fa_file.write(fasta)
    check_run(mod.run_internal_median_filter, fasta, expected_tests[3], str(in_fa_file))

def test_run_internal_median_filter_4(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[4]
    in_fa_file.write(fasta)
    check_run(mod.run_internal_median_filter, fasta, expected_tests[4], str(in_fa_file))

def test_run_internal_median_filter_5(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[5]
    in_fa_file.write(fasta)
    check_run_raises(mod.run_internal_median_filter, fasta_tests[5], str(in_fa_file))

def test_run_internal_median_filter_6(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[6]
    in_fa_file.write(fasta)
    check_run(mod.run_internal_median_filter, fasta_tests[6], expected_tests[6], str(in_fa_file))

def test_run_internal_median_filter_7(tmpdir):
    """
    The run_median_filter asserts if the input file doesn't exist.
    This is to prevent specifying streams we can't rewind.
    """
    check_run_raises(mod.run_internal_median_filter, fasta_tests[5], '-')

def test_run_internal_median_filter_8(tmpdir):
    in_fa_file = tmpdir.join('in.fa')
    fasta = fasta_tests[6]
    in_fa_file.write(fasta)
    check_run(mod.run_internal_median_filter, fasta_tests[7], expected_internal_median_tests[7], str(in_fa_file))


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
