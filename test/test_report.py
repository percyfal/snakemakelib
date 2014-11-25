# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import csv
import unittest
import logging
import texttable as tt
from collections import OrderedDict
from nose.tools import raises
from snakemakelib.report.picard import PicardMetrics, PicardHistMetrics, AlignMetrics, InsertMetrics, HsMetrics, DuplicationMetrics, AlignMetricsCollection, InsertMetricsCollection, HsMetricsCollection, DuplicationMetricsCollection, _read_picard_metrics, PicardMetrics, PicardHist,  PicardMetricsCollection

logger = logging.getLogger(__name__)


class TestCollectMetrics(unittest.TestCase):
    def setUp(self):
        metricsroot = os.path.join(os.path.abspath(os.curdir), os.pardir, 'data', 'metrics', 'J.Doe_00_01')
        self.metricsfiles = []
        for root, dirs, files in os.walk(metricsroot):
            self.metricsfiles += [os.path.join(root, x) for x in files if x.endswith('metrics')]
        self.align_metrics = [x for x in self.metricsfiles if x.endswith('align_metrics')]
        self.hs_metrics = [x for x in self.metricsfiles if x.endswith('hs_metrics')]
        self.insert_metrics = [x for x in self.metricsfiles if x.endswith('insert_metrics')]
        self.dup_metrics = [x for x in self.metricsfiles if x.endswith('dup_metrics')]
        self.alnmet = AlignMetricsCollection(pmid="ID1", file=self.align_metrics[0])
        self.insmet = InsertMetricsCollection(pmid="ID1", file=self.insert_metrics[0])
        self.hsmet = HsMetricsCollection(pmid="ID1", file=self.hs_metrics[0])
        self.dupmet = DuplicationMetricsCollection(pmid="ID1", file=self.dup_metrics[0])

        
    def test_picard_read_metrics(self):
        """Test function _read_picard_metrics"""
        (metrics, hist) = _read_picard_metrics(self.align_metrics[0])
        self.assertIsNone(hist)
        self.assertEqual(len(metrics), 4)
        (metrics, hist) = _read_picard_metrics(self.insert_metrics[0])
        self.assertListEqual(sorted(hist[0]), sorted(['insert_size', 'All_Reads.fr_count']))
        self.assertEqual(metrics[1][0], 156)

    def test_AlignMetrics(self):
        """Test AlignMetrics class"""
        (metrics, hist) = _read_picard_metrics(self.align_metrics[0])
        m = AlignMetrics(*metrics)
        self.assertListEqual(m._fieldnames[0:3], ['CATEGORY', 'TOTAL_READS', 'PF_READS'])
        msum = m.summary(columns=['CATEGORY', 'PCT_PF_READS', 'PF_READS_ALIGNED', 'PCT_PF_READS_ALIGNED'])
        self.assertEqual(msum.split("\n")[1], 'FIRST_OF_PAIR\t1.00\t1.99E+03\t1.00')

    def test_InsertMetrics(self):
        """Test InsertMetrics class"""
        (metrics, hist) = _read_picard_metrics(self.insert_metrics[0])
        i = InsertMetrics(*metrics)
        self.assertListEqual(i.fieldnames[0:3], ['MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE'])
        isum = i.summary(columns=['MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION'])
        self.assertEqual(isum.split("\n")[1], '156\t39')

    def test_HsMetrics(self):
        """Test HsMetrics class"""
        (metrics, hist) = _read_picard_metrics(self.hs_metrics[0])
        h = HsMetrics(*metrics)
        self.assertListEqual(h.fieldnames[0:3],['BAIT_SET', 'GENOME_SIZE', 'BAIT_TERRITORY'])
        hsum = h.summary(columns=['BAIT_SET', 'GENOME_SIZE', 'BAIT_TERRITORY'])
        self.assertEqual(hsum.split("\n")[1], 'chr11_baits\t2.00E+06\t3.01E+02')


    def test_DuplicationMetrics(self):
        """Test DuplicationMetrics class"""
        (metrics, hist) = _read_picard_metrics(self.dup_metrics[0])
        d = DuplicationMetrics(*metrics)
        self.assertListEqual(d.fieldnames[0:3],['LIBRARY', 'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED'])
        dsum = d.summary(columns=['LIBRARY', 'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED'])
        self.assertEqual(dsum.split("\n")[1], 'lib\t5.40E+01\t1.94E+03')

    def test_PicardHist(self):
        """Test PicardHist class"""
        (metrics, hist) = _read_picard_metrics(self.insert_metrics[0])
        h = PicardHist(*hist)
        self.assertListEqual(list(h.hist.keys()), ['insert_size', 'All_Reads.fr_count'])
        self.assertEqual(h.hist['insert_size'][0], 70)
        (metrics, hist) = _read_picard_metrics(self.dup_metrics[0])
        h = PicardHist(*hist)
        self.assertListEqual(list(h.hist.keys()), ['BIN', 'VALUE'])
        self.assertEqual(h.hist['VALUE'][0], 0.99999)

        

    def test_picard_metrics_get_columns(self):
        alnmet = self.alnmet.metrics(columns=['MEAN_READ_LENGTH'])
        print (alnmet)
        self.assertDictEqual(alnmet[0], {'MEAN_READ_LENGTH': '76'})

    @raises(KeyError)
    def test_picard_metrics_get_nonexisting_columns(self):
        alnmet = self.alnmet.metrics(columns=['mMEAN_READ_LENGTH'])
        
    def test_picard_metrics_hist(self):
        insmet = self.insmet.hist()
        self.assertListEqual(sorted(list(insmet.keys())), sorted(['All_Reads.fr_count', 'insert_size']))

    def test_picard_metrics_id(self):
        alnmet = self.alnmet.metrics(columns=['ID', 'MEAN_READ_LENGTH'])
        self.assertDictEqual(alnmet[0], {'ID':'ID1', 'MEAN_READ_LENGTH':'76'})

    # def test_picard_metrics_collection(self):
    #     amc = AlignMetricsCollection([AlignMetrics(pmid=os.path.basename(x).split(".")[0], file=x) for x in self.align_metrics])
    #     am = amc.merge(columns=['ID', 'CATEGORY', 'PCT_PF_READS_ALIGNED'], category=['PAIR'])
    #     self.assertListEqual(am, [{'ID': 'P001_101_index3', 'CATEGORY': 'PAIR', 'PCT_PF_READS_ALIGNED': '0.985015'}, {'ID': 'P001_102_index6', 'CATEGORY': 'PAIR', 'PCT_PF_READS_ALIGNED': '0.981'}])

    def test_align_metrics_category(self):
        """Make sure alignmetrics returns only elements of category"""
        alnmet = self.alnmet.metrics(columns=['ID', 'MEAN_READ_LENGTH', 'CATEGORY'], category=['FIRST_OF_PAIR'])
        self.assertListEqual(alnmet, [{'ID': 'ID1', 'CATEGORY': 'FIRST_OF_PAIR', 'MEAN_READ_LENGTH': '76'}])

    def test_metrics(self):
        print (self.hsmet.metrics())
        print (self.hsmet.hist())
        print (self.dupmet.metrics())
        print (self.dupmet.hist())
        print (self.insmet.metrics())
        print (self.insmet.hist())

    def test_add_metrics(self):
        aml = [AlignMetricsCollection(pmid=os.path.basename(x), file=x) for x in self.align_metrics]
        aml_sum = aml[0] + aml[1]
        print(aml)
        print(aml[0] )
        
    def test_picard_metrics_read_metrics(self):
        """Test PicardMetrics._read_metrics"""
        print (self.alnmet.metrics())
        print (self.insmet.hist())

class TestMergeDict(unittest.TestCase):
    def test_merge_two_dicts(self):
        """Test merging two dictionaries where we know we have the same columns"""
        print ("Merging")

class TestSampleReport(unittest.TestCase):
    def test_sample_report(self):
        """Test sample report"""
        print ("sample report")

