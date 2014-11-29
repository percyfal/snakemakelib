# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import os
import re
import unittest
import logging
import texttable as tt
from collections import OrderedDict, Counter
from nose.tools import raises
from snakemakelib.report.picard import PicardMetrics, PicardHistMetrics, AlignMetrics, InsertMetrics, HsMetrics, DuplicationMetrics, _read_picard_metrics, combine_metrics, DataFrame, _make_unique

logger = logging.getLogger(__name__)

def setUp():
    """Set up test fixtures for metrics test"""

    global PM, PH, AMa, IMa, DMa, HMa, AMf, IMf, DMf, HMf, align_metrics, dup_metrics, insert_metrics, hs_metrics, alnmet, insmet, dupmet, hsmet, inshist, duphist, Adf, Idf, Idfh, Ddf, Ddfh, Hdf

    metricsroot = os.path.join(os.path.abspath(os.curdir),
                               os.pardir, 'data', 'metrics', 'J.Doe_00_01')
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

    Adf = DataFrame(*alnmet)
    Idf = DataFrame(*insmet)
    Idfh = DataFrame(*inshist)
    Ddf = DataFrame(*dupmet)
    Ddfh = DataFrame(*duphist)
    Hdf = DataFrame(*hsmet)

class TestDataFrame(unittest.TestCase):
    """Test DataFrame object"""
    def test_dataframe_init(self):
        """Test initialization of DataFrame"""
        df = DataFrame(*alnmet)
        self.assertListEqual(df.colnames[0:3], ['CATEGORY', 'TOTAL_READS', 'PF_READS'])
        self.assertTupleEqual((3,25), df.dim)
        df = DataFrame(*inshist)
        self.assertTupleEqual((257,2), df.dim)
        self.assertListEqual(df.colnames, ['insert_size', 'All_Reads.fr_count'])

    def test_x(self):
        """Test getting x in two ways"""
        self.assertListEqual(['FIRST_OF_PAIR', 'SECOND_OF_PAIR', 'PAIR'], Adf.x('CATEGORY'))
        self.assertListEqual(['FIRST_OF_PAIR', 'SECOND_OF_PAIR', 'PAIR'], Adf[['CATEGORY']].x())

    def test_y(self):
        """Test getting y in two ways"""
        self.assertListEqual(['FIRST_OF_PAIR', 'SECOND_OF_PAIR', 'PAIR'], Adf.y('CATEGORY'))
        self.assertListEqual(['FIRST_OF_PAIR', 'SECOND_OF_PAIR', 'PAIR'], Adf[['CATEGORY']].y())

    def test_slice_x(self):
        """Test getting x as slice in two ways"""
        self.assertListEqual(Adf.x('CATEGORY', [1]), ['SECOND_OF_PAIR'])
        self.assertListEqual([86, 87, 88, 89, 90, 91, 92, 93, 94, 95], [int(x) for x in Idfh.x(indices=list(range(10,20)))])
        self.assertListEqual([86, 87, 88, 89, 90, 91, 92, 93, 94, 95], [int(x) for x in Idfh.y()[10:20]])
        
    def test_slice_y(self):
        """Test getting y as slice in two ways"""
        self.assertListEqual(Adf.y('CATEGORY', [1]), ['SECOND_OF_PAIR'])
        self.assertListEqual([int(x) for x in Idfh.y(indices=list(range(3,7)))], [79, 80, 81, 82])
        self.assertListEqual([int(x) for x in Idfh.y()[3:7]], [79, 80, 81, 82])

    def test_format(self):
        """Test formatting output of data frame"""
        fmt = OrderedDict([('LIBRARY', ('s', str)), ('UNPAIRED_READS_EXAMINED', ('3.2h', int)), 
                                       ('READ_PAIRS_EXAMINED', ('3.2h', int)), ('UNMAPPED_READS', ('3.2h', int)),
                                       ('UNPAIRED_READ_DUPLICATES', ('3.2h', int)), ('READ_PAIR_DUPLICATES', ('3.2h', int)), 
                                       ('READ_PAIR_OPTICAL_DUPLICATES', ('3.2f', float)), 
                                       ('PERCENT_DUPLICATION', ('3.2%', float)), ('ESTIMATED_LIBRARY_SIZE', ('3.2h', int))])
        df = DataFrame(*dupmet, **fmt)
        self.assertListEqual(sorted(list(df._format.items())), sorted(list(fmt.items())))

class TestReadPicardMetrics(unittest.TestCase):
    """Test reading picard metrics"""

    def test_picard_read_metrics(self):
        """Test function _read_picard_metrics"""
        (metrics, hist) = _read_picard_metrics(align_metrics[0])
        self.assertIsNone(hist)
        self.assertEqual(len(metrics), 4)
        (metrics, hist) = _read_picard_metrics(insert_metrics[0])
        self.assertListEqual(sorted(hist[0]), sorted(['insert_size', 'All_Reads.fr_count']))
        self.assertEqual(metrics[1][0], 156)

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
        self.assertListEqual(p1.metrics.data, p2.metrics.data)

    def test_iter(self):
        """Test PicardMetrics iteration"""
        i = 0
        for m in PM.metrics:
            self.assertListEqual(list(PM.metrics.data[i]), list(m))
            i += 1

    def test_get_single_column(self):
        """Test getting single column"""
        pms = PM.metrics[['SAMPLE']]
        self.assertEqual(pms.colnames, ['SAMPLE'])
        self.assertListEqual([OrderedDict([('SAMPLE', '')]), OrderedDict([('SAMPLE', '')]), OrderedDict([('SAMPLE', '')])], pms.data)

    def test_set_sample(self):
        """Test updating SAMPLE key in PicardMetrics object"""
        PM.metrics['SAMPLE'] = 'sample'
        pms = PM.metrics[['SAMPLE']]
        self.assertListEqual([OrderedDict([('SAMPLE', 'sample')]), OrderedDict([('SAMPLE', 'sample')]), OrderedDict([('SAMPLE', 'sample')])], pms.data)

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
        self.assertListEqual(p1.metrics.data, p2.metrics.data)


class TestAlignMetrics(unittest.TestCase):
    """Test AlignMetrics classes"""
    def test_equality(self):
        """Test that two instances have identical metrics"""
        self.assertListEqual(AMa.metrics.data, AMf.metrics.data)

    def test_summary(self):
        """Test AlignMetrics summary"""
        self.assertEqual('FIRST_OF_PAIR	2.00k	2.00k	100.00%', AMa.summary().split("\n")[1][0:33])

    def test_subset_summary(self):
        """Test AlignMetrics subset summary"""
        columns = ['CATEGORY', 'TOTAL_READS', 'PF_READS_ALIGNED']
        am = AMa[columns]
        self.assertEqual('SECOND_OF_PAIR	2.00k	1.95k', am.summary().split("\n")[2])

    def test_category(self):
        """Test AlignMetrics category retrieval"""
        AMc = AMa.category()
        self.assertTupleEqual(AMc.metrics.dim, (1, 25))
        self.assertEqual(AMc.summary().split("\n")[1].split("\t")[0], "PAIR")
        AMc = AMa.category('FIRST_OF_PAIR')
        self.assertTupleEqual(AMc.metrics.dim, (1, 25))
        self.assertEqual(AMc.summary().split("\n")[1].split("\t")[0], "FIRST_OF_PAIR")


    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass

class TestInsertMetrics(unittest.TestCase):
    """Test InsertMetrics classes"""
    def test_equality(self):
        """Test that two instances have identical metrics"""
        self.assertListEqual(IMa.metrics.data, IMf.metrics.data)

    def test_subset_type(self):
        """Test subsetting an InsertMetrics object"""
        columns = ['MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE']
        im = IMa[columns]
        self.assertListEqual(im.metrics.colnames, columns)
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
        self.assertListEqual(DMa.metrics.data, DMf.metrics.data)

    def test_subset_type(self):
        """Test subsetting an DuplicationMetrics object"""
        columns = ['LIBRARY', 'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED']
        dm = DMa[columns]
        self.assertListEqual(dm.metrics.colnames, columns)
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
        self.assertListEqual(HMa.metrics.data, HMf.metrics.data)

    def test_subset_type(self):
        """Test subsetting an HsMetrics object"""
        columns = ['BAIT_SET', 'GENOME_SIZE', 'BAIT_TERRITORY']
        hm = HMa[columns]
        self.assertListEqual(hm.metrics.colnames, columns)
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


class TestCombineMetrics(unittest.TestCase):
    """Test methods for combining metrics"""

    def test_combine_metrics(self):
        """Test combining metrics"""
        amsa = AMa.category()
        mlist = list(zip([DMa],[IMa],[amsa],[HMa]))
        colnames = [c for sublist in [m.metrics.colnames for mtup in mlist for m in mtup] for c in sublist]
        cnt = Counter(colnames)
        
        # Merge colnames
        pm = combine_metrics(mlist)
        self.assertEqual((len(cnt)), pm.metrics.dim[1])
        # Make colnames unique
        pm = combine_metrics(mlist, uniquify=True)
        self.assertEqual((len(colnames)), pm.metrics.dim[1])

    def test_combine_multiple_metrics(self):
        """Test combining multiple metrics"""
        mlist =(list(
            zip(
                [AlignMetrics(filename=x).category() for x in align_metrics],
                [InsertMetrics(filename=x) for x in insert_metrics],
                [DuplicationMetrics(filename=x) for x in dup_metrics],
                [HsMetrics(filename=x) for x in hs_metrics]
            )
        ))
        pm = combine_metrics(mlist)
        self.assertTupleEqual(pm.metrics.dim, (2,91))
        self.assertListEqual(pm.summary(raw=True).split("\n")[1].split("\t")[:25], [v for k,v in AMa.metrics.data[2].items()])


    def test_plot_tuple(self):
        """Test retrieval of plot tuple"""
        pass
