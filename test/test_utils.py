# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import unittest
import logging
from snakemakelib.config import init_sml_config
from snakemakelib.report.utils import Template
from snakemakelib.bio.ngs.utils import read_group_from_str
from snakemakelib.utils import isoformat

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
        cfg = {
            'bio.ngs.settings' : {
                'run_id_re' : (("platform-unit", "date", "_", "sample"), "([0-9])_([0-9]+)_([A-Z0-9]+XX)_(P[0-9]+_[0-9]+)"), 
                'read_group_keys' : ("id", "sample", "library", "description", "platform-unit", "center", "date", "platform"),
                'center' : 'mycenter',
                'platform' : 'Illumina',
            },
        }
        init_sml_config(cfg)

    def tearDown(self):
        init_sml_config({})

    def test_read_group_from_str(self):
        """Test creating read group information from strings"""
        d = read_group_from_str("2_120924_AC003CCCXX_P001_102")
        self.assertEqual(d['date'], "2012-09-24")
        self.assertEqual(d['platform-unit'], "2")
        self.assertEqual(d['id'], "2_120924_AC003CCCXX_P001_102")

class TestUtils(unittest.TestCase):
    """Test snakemakelib.utils functions"""
    def test_isoformat(self):
        """Test isoformatting function"""
        s = "120924"
        self.assertEqual(isoformat(s), "2012-09-24")
