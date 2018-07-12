"""
This is meant to replace pysiv2.custom.test_assembly.
In your testkit_cfg.json, instead of
    "pysiv2.custom": [
you will now use
    "falcon_kit.testkit": [
The values *might* need to be updated slightly, since we rely on the DB
rather than on datasets now.
"""
import os, re, unittest
from falcon_kit.util.io import system
from falcon_kit.io import capture
from falcon_kit import functional

# Someday, we might simplify this. For now, we
# use the current 'test_values.json'.
from pysiv2.custom.base import TestStatisticsBase # pylint: disable=no-name-in-module, import-error

class TestFalconAssembly(TestStatisticsBase):
    JSON_SCOPE = 'falcon_kit'
    TEST_ID = 'filter_subreads'
    METRIC_IDS = ['number_of_filtered_subreads']
    DEFAULT_VALUES = {}

    @classmethod
    def getMetrics(cls):
        db_fn = os.path.join(cls.job_dir, 'tasks', 'falcon_ns2.tasks.task_falcon0_dazzler_build_raw-0', 'raw_reads.db')
        #system('which DBdump', check=True)
        dump = capture('DBdump {}'.format(db_fn))
        cls.metric_dict['number_of_filtered_subreads'] = functional.dazzler_num_reads(dump)
