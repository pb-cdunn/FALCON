import falcon_kit.mains.fasta_filter as mod
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
>synthetic/2/0_500
A
>synthetic/3/0_500
A
>synthetic/1/0_500
ACGTA
>synthetic/1/0_500
AC
>synthetic/1/0_500
ACG
>synthetic/1/0_500
ACGT
>synthetic/1/0_500
A
>synthetic/4/0_500
A
>synthetic/5/0_500
AC
>synthetic/5/0_500
ACG
>synthetic/5/0_500
ACGT
>synthetic/5/0_500
ACGTA
"""
expected_tests[4] = """\
>synthetic/2/0_500
A
>synthetic/3/0_500
A
>synthetic/1/0_500
ACG
>synthetic/4/0_500
A
>synthetic/5/0_500
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
>synthetic/5/0_500
ACG
>synthetic/2/0_500
A
>synthetic/3/0_500
A
>synthetic/1/0_500
ACGTA
>synthetic/1/0_500
ACGT
>synthetic/1/0_500
A
>synthetic/4/0_500
A
>synthetic/5/0_500
AC
>synthetic/5/0_500
ACGT
>synthetic/5/0_500
ACGTA
>synthetic/1/0_500
AC
>synthetic/1/0_500
ACG
"""
expected_tests[6] = """\
>synthetic/2/0_500
A
>synthetic/3/0_500
A
>synthetic/4/0_500
A
>synthetic/5/0_500
ACGT
>synthetic/1/0_500
ACG
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
    check_run(mod.run_streamed_median, fasta_tests[1], expected_tests[1], '-')

def test_run_streamed_median_2():
    check_run(mod.run_streamed_median, fasta_tests[2], expected_tests[2], '-')

def test_run_streamed_median_3():
    check_run(mod.run_streamed_median, fasta_tests[3], expected_tests[3], '-')

def test_run_streamed_median_4():
    check_run(mod.run_streamed_median, fasta_tests[4], expected_tests[4], '-')

def test_run_streamed_median_5():
    check_run_raises(mod.run_streamed_median, fasta_tests[5], '-')

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

def test_run_median_filter_7(tmpdir):
    """
    The run_median_filter asserts if the input file doesn't exist.
    This is to prevent specifying streams we can't rewind.
    """
    check_run_raises(mod.run_median_filter, fasta_tests[5], '-')

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
