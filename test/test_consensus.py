from cStringIO import StringIO
import sys
import falcon_kit.mains.consensus as mod

class Dummy(object): pass

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_process_get_consensus_result(capsys):
    res = ('A', 0)
    args = Dummy()
    mod.process_get_consensus_result(res, args)
    (out, err) = capsys.readouterr()
    assert out == ''

    res = ('abcd', '42')
    args = Dummy()
    args.output_full = True
    mod.process_get_consensus_result(res, args, limit=4)
    (out, err) = capsys.readouterr()
    assert out == '>42_f\nabcd\n'

    res = ('nope', 0)
    args = Dummy()
    args.output_full = False
    args.output_multi = False
    mod.process_get_consensus_result(res, args, limit=3)
    (out, err) = capsys.readouterr()
    assert out == ''

    res = ('nope', 0)
    args = Dummy()
    args.output_full = False
    args.output_multi = True
    mod.process_get_consensus_result(res, args, limit=3)
    (out, err) = capsys.readouterr()
    assert out == ''

    res = ('ATATXCGC', 77)
    args = Dummy()
    args.output_full = False
    args.output_multi = False
    mod.process_get_consensus_result(res, args, limit=2)
    (out, err) = capsys.readouterr()
    assert out == '>77\nATAT\n'

    base = 'AA'
    full = 'x'.join(['G'] + [base]*11 + ['T'])
    res = (full, 1)
    args = Dummy()
    args.output_full = False
    args.output_multi = True
    mod.process_get_consensus_result(res, args, limit=2)
    (out, err) = capsys.readouterr()
    lines = out.split('\n')
    expected = [
        u'>prolog/10/0_2',
        u'AA',
        u'>prolog/11/0_2',
        u'AA',
        u'>prolog/12/0_2',
        u'AA',
        u'>prolog/13/0_2',
        u'AA',
        u'>prolog/14/0_2',
        u'AA',
        u'>prolog/15/0_2',
        u'AA',
        u'>prolog/16/0_2',
        u'AA',
        u'>prolog/17/0_2',
        u'AA',
        u'>prolog/18/0_2',
        u'AA',
        u'>prolog/19/0_2',
        u'AA',
        u'',
    ]
    assert lines == expected
