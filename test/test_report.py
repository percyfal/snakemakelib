# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import os
import unittest
import logging
from collections import namedtuple
from snakemake.utils import report
from pylab import *
import matplotlib.pyplot as plt

from snakemakelib.report.picard import PicardMetrics, PicardHistMetrics, AlignMetrics, InsertMetrics, HsMetrics, DuplicationMetrics, combine_metrics

logger = logging.getLogger(__name__)

def setUp():
    """Set up test fixtures"""

    global mlist, Sample, pm

    Sample = namedtuple('Sample', ['sample_id', 'project_id', 'pm'])
    
    metricsroot = os.path.join(os.path.abspath(os.curdir),
                               os.pardir, 'data', 'metrics', 'J.Doe_00_01')
    metricsfiles = []
    for root, dirs, files in os.walk(metricsroot):
        metricsfiles += [os.path.join(root, x) for x in files if x.endswith('metrics')]

    align_metrics = [x for x in metricsfiles if x.endswith('align_metrics')]
    hs_metrics = [x for x in metricsfiles if x.endswith('hs_metrics')]
    insert_metrics = [x for x in metricsfiles if x.endswith('insert_metrics')]
    dup_metrics = [x for x in metricsfiles if x.endswith('dup_metrics')]

    mlist =(list(
        zip(
            [AlignMetrics(filename=x).category() for x in align_metrics],
            [InsertMetrics(filename=x) for x in insert_metrics],
            [DuplicationMetrics(filename=x) for x in dup_metrics],
            [HsMetrics(filename=x) for x in hs_metrics]
        )
    ))
    #pm = combine_metrics(mlist)

class TestSampleReport(unittest.TestCase):
    """Test sample reports"""
    def test_index(self):
        """Test the index page"""
        pass

    def test_sample_report(self):
        """Test sample report"""
        pass
