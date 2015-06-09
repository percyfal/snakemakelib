# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import io
import pandas as pd
from unittest.mock import patch
from snakemakelib.bio.ngs.align.star import Star


class TestStar(unittest.TestCase):
    """Test star"""
    def setUp(self):
        self.data = """                                 Started job on |       Apr 29 22:42:00
                             Started mapping on |       Apr 29 22:42:00
                                    Finished on |       Apr 29 22:43:51
       Mapping speed, Million of reads per hour |       152.43

                          Number of input reads |       4699845
                      Average input read length |       202
                                    UNIQUE READS:
                   Uniquely mapped reads number |       4011114
                        Uniquely mapped reads % |       85.35%
                          Average mapped length |       198.26
                       Number of splices: Total |       1452777
            Number of splices: Annotated (sjdb) |       1424534
                       Number of splices: GT/AG |       1429760
                       Number of splices: GC/AG |       12299
                       Number of splices: AT/AC |       1528
               Number of splices: Non-canonical |       9190
                      Mismatch rate per base, % |       0.70%
                         Deletion rate per base |       0.02%
                        Deletion average length |       1.76
                        Insertion rate per base |       0.01%
                       Insertion average length |       1.46
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       267393
             % of reads mapped to multiple loci |       5.69%
        Number of reads mapped to too many loci |       7530
             % of reads mapped to too many loci |       0.16%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       8.73%
                     % of reads unmapped: other |       0.08%
"""
        self.f = io.StringIO(self.data)

    def test_collect_results(self):
        st = Star([(self.f, 'bar')])
        self.assertEqual(st['align'].shape[0], 1)
        self.assertEqual(st['align'].shape[1], 27)
