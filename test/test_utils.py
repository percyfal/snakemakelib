# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import unittest
import logging
from nose.tools import raises
from snakemakelib.config import init_sml_config
from snakemakelib.report.utils import Template
from snakemakelib.bio.ngs.utils import find_files, ReadGroup, DisallowedKeyException
from snakemakelib.utils import isoformat
from snakemake.io import glob_wildcards

logger = logging.getLogger(__name__)

class TestReportUtils(unittest.TestCase):
    """Test report utilities"""
    def test_format(self):
        """Test number formatting function"""
        tp = Template()
        T = 12121414141232
        G = 12121414141.235
        M = 12121414
        k = 12123.23
        z = 12.25
        m = 0.01212
        u = 0.00001212
        n = 0.00000001212
        Ts = tp.format_field(T, '3.2h')
        Gs = tp.format_field(G, '3.2h')
        Ms = tp.format_field(M, '3.2h')
        ks = tp.format_field(k, '3.2h')
        zs = tp.format_field(z, '3.2h')
        ms = tp.format_field(m, '3.2h')
        us = tp.format_field(u, '3.2h')
        ns = tp.format_field(n, '3.2h')
        self.assertEqual(Ts, '12.12T')
        self.assertEqual(Gs, '12.12G')
        self.assertEqual(Ms, '12.12M')
        self.assertEqual(ks, '12.12k')
        self.assertEqual(zs, '12.25')
        self.assertEqual(ms, '12.12m')
        self.assertEqual(us, '12.12u')
        self.assertEqual(ns, '12.12n')

class TestBioNgsUtils(unittest.TestCase):
    """Test snakemakelib.bio.ngs.utils functions"""
    def setUp(self):
        """Setup test fixtures"""
        self.cfg = {
            'bio.ngs.settings' : {
                'run_id_re' : (("punit", "date", "_", "sample"), "([0-9])_([0-9]+)_([A-Z0-9]+XX)_(P[0-9]+_[0-9]+)"), 
                'read_group_keys' : ("id", "sample", "library", "description", "punit", "center", "date", "platform"),
                'center' : 'mycenter',
                'platform' : 'Illumina',
                'fastq_suffix' : ".fastq.gz",
                'read1_label' : "_1",
                'read2_label' : "_2",
            },
        }
        init_sml_config(self.cfg)

    def tearDown(self):
        init_sml_config({})

    def test_find_files(self):
        """Test finding files"""
        read1_str = self.cfg['bio.ngs.settings']["run_id_re"][1] + self.cfg['bio.ngs.settings']["read1_label"] + self.cfg['bio.ngs.settings']["fastq_suffix"] + "$"
        read2_str = self.cfg['bio.ngs.settings']["run_id_re"][1] + self.cfg['bio.ngs.settings']["read2_label"] + self.cfg['bio.ngs.settings']["fastq_suffix"] + "$"
        flist = find_files("../data/projects/J.Doe_00_01", read1_str)
        self.assertEqual(len(flist), 3)
        self.assertEqual(flist[1], '../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz')
        readpairs = list(zip(flist, find_files("../data/projects/J.Doe_00_01", read2_str)))
        self.assertEqual(len(readpairs), 3)
        self.assertEqual(len(readpairs[0]), 2)
        self.assertEqual(readpairs[1][0], '../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz')

    def test_sm_wildcard(self):
        read1_str = "../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/" + self.cfg['bio.ngs.settings']["run_id_re"][1] + self.cfg['bio.ngs.settings']["read1_label"] + self.cfg['bio.ngs.settings']["fastq_suffix"] + "$"
        print (read1_str)
        gw = glob_wildcards(read1_str)
        print (gw)
        read2_str = "../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/{platform_unit}_{fc}_{sample}_2.fastq.gz"
        gw = glob_wildcards(read2_str)
        print (gw)
        

class TestUtils(unittest.TestCase):
    """Test snakemakelib.utils functions"""
    def test_isoformat(self):
        """Test isoformatting function"""
        s = "120924"
        self.assertEqual(isoformat(s), "2012-09-24")

class TestReadGroup(unittest.TestCase):
    """Test ReadGroup class"""
    def test_rg_init(self):
        """Test initializing ReadGroup"""
        rg = ReadGroup("test", ID='test', DT="120924")
        self.assertEqual(str(rg), '--date 2012-09-24 --identifier test')
        rg = ReadGroup("test", **{'ID':'test', 'DT':"120924"})
        self.assertEqual(str(rg), '--date 2012-09-24 --identifier test')

    def test_rg_parse_illumina_like(self):
        """Test parsing illumina-like-based file names"""
        rg = ReadGroup("(?P<PATH>.*)(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)")
        s = (rg.parse("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz"))
        self.assertEqual("--date 2012-10-15 --identifier 1_121015_BB002BBBXX_P001_101 --platform-unit 1_BB002BBBXX --sample P001_101", s)

    def test_rg_fn(self):
        """Test initializing read group class and setting function"""
        rg_fn = ReadGroup("(?P<PATH>.*)(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)").parse
        s = rg_fn("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz")
        self.assertEqual("--date 2012-10-15 --identifier 1_121015_BB002BBBXX_P001_101 --platform-unit 1_BB002BBBXX --sample P001_101", s)

    @raises(DisallowedKeyException)
    def test_rg_disallowed_key(self):
        """Test setting a read group object with a key not present in allowed keys"""
        rg = ReadGroup("(?P<PATH>.*)(?P<PU1>[0-9])_(?P<DATE>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)")
        s = (rg.parse("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz"))
