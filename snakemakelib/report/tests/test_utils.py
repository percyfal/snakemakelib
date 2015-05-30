# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import datetime
from snakemakelib.report.utils import Template, recast

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

    def test_recast(self):
        """Test recasting string to float, int, date, or other"""
        self.assertEqual(type(recast("1234")), int)
        self.assertEqual(type(recast("123.45")), float)
        self.assertEqual(type(recast("123,45")), float)
        self.assertEqual(type(recast("23.45%")), float)
        self.assertEqual(type(recast("23,45%")), float)
        self.assertEqual(type(recast("Mar 23 00:24:12")), datetime.datetime)

