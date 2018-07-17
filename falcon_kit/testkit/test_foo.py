import os, unittest
from falcon_kit.util.io import system

# Also in pbcommand/testkit/base_utils.py
def pb_requirements(*reqs):
    """
    Method decorator for specifying linked JIRA issues.
    """
    def decorator(test_item):
        test_item.__pb_requirements__ = list(reqs)
        return test_item
    return decorator

class TestBase(unittest.TestCase):
    _is_test = True
    job_dir = None
    service_access_layer = None
    job_id = None

@unittest.skip('FOO')
class TestMe(TestBase):
    @unittest.skip('would fail')
    def test_excfunc(self):
        raise Exception('FAIL HERE')
    @unittest.skip('would fail')
    def test_assertfunc(self):
        assert 0 > 1, 'no way'
    @unittest.skip('would fail')
    def test_failfunc(self):
        msg = 'In {} dir'.format(self.job_dir)
        self.fail(msg)
    def test_goodfunc(self):
        pass
    @pb_requirements('TAGT-000')
    def test_myprd_func(self):
        pass
