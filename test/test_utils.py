# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import unittest
import math
import logging
from nose.tools import raises
from snakemakelib.report.utils import Template

logger = logging.getLogger(__name__)

class TestUtils(unittest.TestCase):
    def test_format(self):
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


        
