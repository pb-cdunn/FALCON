from falcon_kit import FastaReader as M
from falcon_kit.io import NativeIO as StringIO

def test_fasta_empty():
    fasta = ''
    sin = StringIO(fasta)
    result = list(M.yield_fasta_record(sin, None))
    assert not result

def test_fasta_simple():
    fasta = """\
>foo/bar/0_42 FOO=BAR
ACGT
"""
    sin = StringIO(fasta)
    result = list(M.yield_fasta_record(sin, None))
    assert [rec.sequence for rec in result] == ['ACGT']
    rec = result[0]
    assert rec.metadata == 'FOO=BAR'
    assert rec.id == 'foo/bar/0_42'
    assert rec.name == 'foo/bar/0_42 FOO=BAR'
