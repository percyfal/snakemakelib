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

    AMa = AlignMetrics(*alnmet)
    IMa = InsertMetrics(*insmet, hist=inshist)
    HMa = HsMetrics(*hsmet)
    DMa = DuplicationMetrics(*dupmet, hist=duphist)


class TestPicardMetrics(unittest.TestCase):
    """Test PicardMetrics classes"""

    @raises(ValueError)
    def test_missing_args(self):
        """Test instantiating PicardMetrics with missing arguments"""
        PicardMetrics()

    def test_init(self):
        """Test instantiating PicardMetrics in two ways"""
        p1 = PicardMetrics(filename=align_metrics[0])
        p2 = PicardMetrics(*alnmet)
        self.assertListEqual(p1.metrics, p2.metrics)

    def test_iter(self):
        """Test PicardMetrics iteration"""
        i = 0
        for m in PM:
            self.assertListEqual(list(PM.metrics[i]), list(m))
            i += 1


    def test_PicardMetrics_merge(self):
        """Test merging two PicardMetrics objects"""
        pass

    def test_PicardMetrics_merge_subset(self):
        """Test merging two subsetted PicardMetrics objects"""
        pass

class TestPicardHistMetrics(unittest.TestCase):
    
    @raises(ValueError)
    def test_PicardHistMetrics_missing_args(self):
        """Test instantiating PicardHistMetrics with missing arguments"""
        PicardHistMetrics()

    @raises(ValueError)
    def test_PicardHistMetrics_missing_hist(self):
        """Test instantiating PicardHistMetrics with missing hist argument"""
        args = [('MEDIAN_INSERT_SIZE', '156'), ('MEDIAN_ABSOLUTE_DEVIATION', '39')]
        p = PicardHistMetrics(*args)

    def test_PicardHistMetrics_init(self):
        """Test instantiating PicardHistMetrics object"""
        p1 = PicardHistMetrics(filename=insert_metrics[0], hist="test")
        p2 = PicardHistMetrics(*insmet, hist="test")
        self.assertListEqual(p1.metrics, p2.metrics)


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
        print (AMa.summary())

    def test_subset_summary(self):
        """Test AlignMetrics subset summary"""
        columns = ['CATEGORY', 'TOTAL_READS', 'PF_READS']
        am = AMa[columns]
        print (am.summary())

    def test_merge(self):
        """Test merging two AlignMetrics objects"""
        pass

    def test_merge_subset(self):
        """Test merging two subsetted AlignMetrics objects"""
        pass

    # @raises
    def test_merge_subset_different_columns(self):
        """Test merging two subsetted AlignMetrics objects with different columns"""
        pass

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
        print (IMa.summary())

    def test_subset_summary(self):
        """Test InsertMetrics subset summary"""
        columns = ['MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE']
        im = IMa[columns]
        print (im.summary())

    def test_merge(self):
        """Test merging two InsertMetrics objects"""
        pass

    def test_merge_subset(self):
        """Test merging two subsetted InsertMetrics objects"""
        pass

    # @raises
    def test_merge_subset_different_columns(self):
        """Test merging two subsetted InsertMetrics objects with different columns"""
        pass

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
        print (DMa.summary())

    def test_subset_summary(self):
        """Test DuplicationMetrics subset summary"""
        columns = ['LIBRARY', 'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED']
        dm = DMa[columns]
        print (dm.summary())

    def test_merge(self):
        """Test merging two DuplicationMetrics objects"""
        pass

    def test_merge_subset(self):
        """Test merging two subsetted DuplicationMetrics objects"""
        pass

    # @raises
    def test_merge_subset_different_columns(self):
        """Test merging two subsetted DuplicationMetrics objects with different columns"""
        pass

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
        print (HMa.summary())

    def test_subset_summary(self):
        """Test HsMetrics subset summary"""
        columns = ['BAIT_SET', 'GENOME_SIZE', 'BAIT_TERRITORY']
        hm = HMa[columns]
        print (hm.summary())

    def test_merge(self):
        """Test merging two HsMetrics objects"""
        pass

    def test_merge_subset(self):
        """Test merging two subsetted HsMetrics objects"""
        pass

    # @raises
    def test_merge_subset_different_columns(self):
        """Test merging two subsetted HsMetrics objects with different columns"""
        pass

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

