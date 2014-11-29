# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import os
import unittest
import logging
from nose.tools import raises
from collections import namedtuple
from mako.template import Template
from snakemake.utils import report
import matplotlib
#from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
import matplotlib.pyplot as plt

from snakemakelib.report.picard import PicardMetrics, PicardHistMetrics, AlignMetrics, InsertMetrics, HsMetrics, DuplicationMetrics, combine_metrics

logger = logging.getLogger(__name__)

TEMPLATEPATH = os.path.join(os.path.dirname(__file__), os.pardir, "snakemakelib", "data", "templates", "doc")

TEMPLATES = {
    'make' : Template(filename=os.path.join(TEMPLATEPATH, "Makefile.mako")),
    'sample' : Template(filename=os.path.join(TEMPLATEPATH, "source", "samples", "sample.mako")),
    'index' : Template(filename=os.path.join(TEMPLATEPATH, "source", "index.mako")),
    'sampleindex' : Template(filename=os.path.join(TEMPLATEPATH, "source", "samples", "index.mako")),
    'conf' : Template(filename=os.path.join(TEMPLATEPATH, "source", "conf.mako")),
}

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
    pm = combine_metrics(mlist)

class TestSampleReport(unittest.TestCase):
    """Test sample reports"""
    def test_index(self):
        """Test the index page"""
        pct_aligned = ([100 * float(x) for x in pm.metrics.x('PCT_PF_READS_ALIGNED')])
        nseq = ([int(x) for x in pm.metrics.x('TOTAL_READS')])
        sdup = ([100 * float(x) for x in pm.metrics.x('PERCENT_DUPLICATION')])
        plt.scatter(pct_aligned, nseq, s=sdup, alpha=0.75)
        plt.xlabel(r'Percent aligned', fontsize=14)
        plt.yscale('log', **{'basey':10})
        plt.ylabel(r'Read count', fontsize=14)
        plt.title("Sequence summary.\nPoint sizes correspond to duplication levels.", fontsize=14)
        plt.tight_layout()
        savefig('test.png')


    def test_sample_report(self):
        """Test sample report"""
        pct_aligned = ([100 * float(x) for x in pm.metrics.x('PCT_PF_READS_ALIGNED')])
        nseq = ([int(x) for x in pm.metrics.x('TOTAL_READS')])
        sdup = ([100 * float(x) for x in pm.metrics.x('PERCENT_DUPLICATION')])
        input = {'test' : 'test.png'}
        kw = {'pct_aligned' : pct_aligned,
              'nseq' : nseq,
              'sdup' : sdup,
              'project_name' : 'J.Doe_00_01',
              'application' : 'application',
              'date' : 'date',
        }
        kw.update(input)

        tmp = TEMPLATES['index'].render(**kw)
        print (tmp)

        report("""
Project summary
=============================

:Project: {project_name}
:Application: {application}
:Date: {date}

Samples
--------

QC Metrics
----------

Sequence statistics
^^^^^^^^^^^^^^^^^^^

.. figure:: {test}



        """.format(**kw), "test.html", **input)


