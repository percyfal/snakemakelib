# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import csv
import unittest
import logging
import texttable as tt
from collections import OrderedDict
from nose.tools import raises
from snakemakelib.report.picard import PicardMetrics, PicardHistMetrics, AlignMetrics, InsertMetrics, HsMetrics, DuplicationMetrics, _read_picard_metrics

logger = logging.getLogger(__name__)

def setUp():
    """Set up test fixtures for metrics test"""

    global PM, PH, AMa, IMa, DMa, HMa, AMf, IMf, DMf, HMf, align_metrics, dup_metrics, insert_metrics, hs_metrics, alnmet, insmet, dupmet, hsmet, inshist, duphist

    metricsroot = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'metrics', 'J.Doe_00_01')
    metricsfiles = []
    for root, dirs, files in os.walk(metricsroot):
        metricsfiles += [os.path.join(root, x) for x in files if x.endswith('metrics')]

    align_metrics = [x for x in metricsfiles if x.endswith('align_metrics')]
    hs_metrics = [x for x in metricsfiles if x.endswith('hs_metrics')]
    insert_metrics = [x for x in metricsfiles if x.endswith('insert_metrics')]
    dup_metrics = [x for x in metricsfiles if x.endswith('dup_metrics')]

    (alnmet, _) = _read_picard_metrics(align_metrics[0])
    (insmet, inshist) = _read_picard_metrics(insert_metrics[0])
    (hsmet, _) = _read_picard_metrics(hs_metrics[0])
    (dupmet, duphist) = _read_picard_metrics(dup_metrics[0])

    PM = PicardMetrics(identifier="PM", filename=align_metrics[0])
    PH = PicardHistMetrics(identifier="PH", filename=insert_metrics[0])
    AMf = AlignMetrics(identifier="AM", filename=align_metrics[0])
    IMf = InsertMetrics(identifier="IM", filename=insert_metrics[0])
    HMf = HsMetrics(identifier="HM", filename=hs_metrics[0])
    DMf = DuplicationMetrics(identifier="DM", filename=dup_metrics[0])

    AMa = AlignMetrics(*alnmet, identifier="AM")
    IMa = InsertMetrics(*insmet, hist=inshist, identifier="IM")
    HMa = HsMetrics(*hsmet, identifier="HM")
    DMa = DuplicationMetrics(*dupmet, hist=duphist, identifier="DM")


class TestPicardMetrics(unittest.TestCase):
    """Test PicardMetrics classes"""

    @raises(ValueError)
    def test_missing_args(self):
        """Test instantiating PicardMetrics with missing arguments"""
        PicardMetrics()

    def test_init(self):
        """Test instantiating PicardMetrics in two ways"""
        p1 = PicardMetrics(filename=align_metrics[0], identifier="PM")
        p2 = PicardMetrics(*alnmet, identifier="PM")
        self.assertListEqual(p1.metrics, p2.metrics)

    def test_iter(self):
        """Test PicardMetrics iteration"""
        i = 0
        for m in PM:
            self.assertListEqual(list(PM.metrics[i]), list(m))
            i += 1

    def test_merge(self):
        """Test merging two PicardMetrics objects"""
        PM2 = PicardMetrics(identifier="PM2", filename=align_metrics[1])
        self.assertEqual(len(list(PM.metrics)), 3)
        self.assertEqual(len(list(PM2.metrics)), 3)
        PM3 = PM + PM2
        self.assertEqual(len(list(PM3.metrics)), 6)

        pm1_cat = [row['CATEGORY'] for row in PM.metrics]
        pm2_cat = [row['CATEGORY'] for row in PM2.metrics]
        pm3_cat = [row['CATEGORY'] for row in PM3.metrics]
        
        pm1_pct = [row['PCT_READS_ALIGNED_IN_PAIRS'] for row in PM.metrics]
        pm2_pct = [row['PCT_READS_ALIGNED_IN_PAIRS'] for row in PM2.metrics]
        pm3_pct = [row['PCT_READS_ALIGNED_IN_PAIRS'] for row in PM3.metrics]
        self.assertListEqual(pm3_cat, pm1_cat + pm2_cat)
        self.assertListEqual(pm3_pct, pm1_pct + pm2_pct)

    def test_merge_subset(self):
        """Test merging two subsetted PicardMetrics objects"""
        PM2 = PicardMetrics(identifier="PM2", filename=align_metrics[1])
        columns = ['CATEGORY', 'PCT_PF_READS_ALIGNED']
        p1 = PM[columns]
        p2 = PM2[columns]
        self.assertEqual(len(list(p1.metrics)), 3)
        self.assertEqual(len(list(p2.metrics)), 3)
        p3 = p1 + p2
        self.assertEqual(len(list(p3.metrics)), 6)
        self.assertListEqual(p3.metrics, p1.metrics + p2.metrics)

    @raises(TypeError)
    def test_merge_different_columns(self):
        """Test merging two PicardMetrics objects with different columns"""
        (pm1, pm2) = (PicardMetrics(identifier="PM1", filename=align_metrics[0]), PicardMetrics(identifier="PM2", filename=align_metrics[1]))
        pm1s = pm1[['CATEGORY']]
        pm1s + pm2

    def test_get_single_column(self):
        """Test getting single column"""
        pms = PM[['SAMPLE']]
        self.assertEqual(pms.fieldnames, ['SAMPLE'])
        self.assertListEqual([OrderedDict([('SAMPLE', '')]), OrderedDict([('SAMPLE', '')]), OrderedDict([('SAMPLE', '')])], pms.metrics)

    def test_set_sample(self):
        """Test updating SAMPLE key in PicardMetrics object"""
        PM['SAMPLE'] = 'sample'
        pms = PM[['SAMPLE']]
        self.assertListEqual([OrderedDict([('SAMPLE', 'sample')]), OrderedDict([('SAMPLE', 'sample')]), OrderedDict([('SAMPLE', 'sample')])], pms.metrics)

class TestPicardHistMetrics(unittest.TestCase):
    
    @raises(ValueError)
    def test_missing_args(self):
        """Test instantiating PicardHistMetrics with missing arguments"""
        PicardHistMetrics()

    @raises(ValueError)
    def test_missing_hist(self):
        """Test instantiating PicardHistMetrics with missing hist argument"""
        args = [('MEDIAN_INSERT_SIZE', '156'), ('MEDIAN_ABSOLUTE_DEVIATION', '39')]
        p = PicardHistMetrics(*args)

    def test_init(self):
        """Test instantiating PicardHistMetrics object"""
        p1 = PicardHistMetrics(filename=insert_metrics[0], hist="test", identifier="PM")
        p2 = PicardHistMetrics(*insmet, hist="test", identifier="PM")
        self.assertListEqual(p1.metrics, p2.metrics)

    def test_merge(self):
        """Test merging two PicardHistMetrics objects"""
        pass

    def test_merge_subset(self):
        """Test merging two subsetted PicardHistMetrics objects"""
        pass

    # @raises
    def test_merge_subset_different_columns(self):
        """Test merging two subsetted PicardHistMetrics objects with different columns"""
        pass


class TestAlignMetrics(unittest.TestCase):
    """Test AlignMetrics classes"""
    def test_equality(self):
        """Test that two instances have identical metrics"""
        self.assertListEqual(AMa.metrics, AMf.metrics)

    def test_subset_type(self):
        """Test subsetting an AlignMetrics object"""
        columns = ['CATEGORY', 'TOTAL_READS', 'PF_READS']
        am = AMa[columns]
        self.assertListEqual(am.fieldnames, columns)
        self.assertIsInstance(am, AlignMetrics)

    def test_summary(self):
        """Test AlignMetrics summary"""
        self.assertEqual('FIRST_OF_PAIR	2.00k	2.00k	100.00%', AMa.summary().split("\n")[1][0:33])

    def test_subset_summary(self):
        """Test AlignMetrics subset summary"""
        columns = ['CATEGORY', 'TOTAL_READS', 'PF_READS_ALIGNED']
        am = AMa[columns]
        self.assertEqual('SECOND_OF_PAIR	2.00k	1.95k', am.summary().split("\n")[2])

    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass

class TestInsertMetrics(unittest.TestCase):
    """Test InsertMetrics classes"""
    def test_equality(self):
        """Test that two instances have identical metrics"""
        self.assertListEqual(IMa.metrics, IMf.metrics)

    def test_subset_type(self):
        """Test subsetting an InsertMetrics object"""
        columns = ['MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE']
        im = IMa[columns]
        self.assertListEqual(im.fieldnames, columns)
        self.assertIsInstance(im, InsertMetrics)

    def test_summary(self):
        """Test InsertMetrics summary"""
        self.assertListEqual(IMa.summary().split("\n")[1].split("\t"), ['156', '39', '70', '485', '167.819', '61.549', '1.73k', 'FR', '15', '29', '43', '61', '79', '93', '111', '133', '195', '443', '', '', ''])

    def test_subset_summary(self):
        """Test InsertMetrics subset summary"""
        columns = ['MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE']
        im = IMa[columns]
        self.assertListEqual(im.summary().split("\n")[1].split("\t"), ['156', '39', '70'])

    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass


class TestDuplicationMetrics(unittest.TestCase):
    """Test DuplicationMetrics classes"""
    def test_equality(self):
        """Test that two instances have identical metrics"""
        self.assertListEqual(DMa.metrics, DMf.metrics)

    def test_subset_type(self):
        """Test subsetting an DuplicationMetrics object"""
        columns = ['LIBRARY', 'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED']
        dm = DMa[columns]
        self.assertListEqual(dm.fieldnames, columns)
        self.assertIsInstance(dm, DuplicationMetrics)

    def test_summary(self):
        """Test DuplicationMetrics summary"""
        self.assertEqual(DMa.summary().split("\n")[1].split("\t"), ['lib', '54.00', '1.95k', '60.00', '43.00', '215.00', '0.00', '11.99%', '8.14k'])

    def test_subset_summary(self):
        """Test DuplicationMetrics subset summary"""
        columns = ['LIBRARY', 'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED', 'PERCENT_DUPLICATION']
        dm = DMa[columns]
        self.assertListEqual(['lib', '54.00', '1.95k', '11.99%'], dm.summary().split("\n")[1].split("\t"))

    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass


class TestHsMetrics(unittest.TestCase):
    """Test HsMetrics classes"""
    def test_equality(self):
        """Test that two instances have identical metrics"""
        self.assertListEqual(HMa.metrics, HMf.metrics)

    def test_subset_type(self):
        """Test subsetting an HsMetrics object"""
        columns = ['BAIT_SET', 'GENOME_SIZE', 'BAIT_TERRITORY']
        hm = HMa[columns]
        self.assertListEqual(hm.fieldnames, columns)
        self.assertIsInstance(hm, HsMetrics)

    def test_summary(self):
        """Test HsMetrics summary"""
        self.assertListEqual(HMa.summary().split("\n")[1].split("\t")[0:6], ['chr11_baits', '2.00M', '301.00', '301.00', '1.00', '4.00k'])

    def test_subset_summary(self):
        """Test HsMetrics subset summary"""
        columns = ['GENOME_SIZE', 'BAIT_TERRITORY', 'TARGET_TERRITORY', 'ZERO_CVG_TARGETS_PCT', 'PCT_TARGET_BASES_2X', 'PCT_TARGET_BASES_10X', 'PCT_TARGET_BASES_30X']
        hm = HMa[columns]
        self.assertListEqual(['2.00M', '301.00', '301.00', '0.00%', '76.41%', '50.83%', '16.61%'], hm.summary().split("\n")[1].split("\t"))

    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass


class TestPicardMetricsSummary(unittest.TestCase):
    """Test PicardMetricsSummary class"""
    def test_summary(self):
        """Test summary combining all metrics, all AlignMetrics rows"""
        pass

    def test_summary_pairs(self):
        """Test summary using only PAIRS category of AlignMetrics"""
        pass

    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass

