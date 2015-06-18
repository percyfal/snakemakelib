# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import pandas as pd
from unittest.mock import patch
from snakemakelib.bio.ngs.qc.qualimap import Qualimap


class TestQualimap(unittest.TestCase):
    """Test Qualimap"""
    def setUp(self):
        self.data = ['>>>>>>> Globals\n',
                     'number of windows = 10\n',
                     'number of reads = 10,000,000\n',
                     'number of mapped reads = 9,900,000 (99.00%)\n',
                     'number of duplicated reads = 400,000\n\n',
                     '>>>>>>> Insert size\n',
                     '>>>>>>> Coverage per contig\n',
                     '\n',
                     '\t'.join(['foo', '11', '12', '1.1', '1.2\n']),
                     '\t'.join(['bar', '21', '22', '2.1', '2.2\n'])]

    @patch('snakemakelib.results.Results.load_lines')
    def test_qualimap_coverage(self, mock_load_lines):
        mock_load_lines.return_value = self.data
        qm = Qualimap([('foo', 'bar')])
        self.assertListEqual([34.375, 65.625],
                             list(qm['coverage_per_contig']['chrlen_percent']))

    @patch('snakemakelib.results.Results.load_lines')
    def test_qualimap_globals(self, mock_load_lines):
        mock_load_lines.return_value = self.data
        qm = Qualimap([('foo', 'bar')])
        self.assertListEqual(sorted([10.0, 400000.0, 9500000.0, 10000000.0, 9900000.0]),
                             sorted(list(qm['globals'].loc['bar'])))

    @patch('pandas.DataFrame')
    @patch('snakemakelib.results.Results.load_lines')
    def test_collect_results(self, mock_load_lines, mock_df):
        mock_load_lines.return_value = self.data
        mock_df.return_value = pd.DataFrame()
        Qualimap([('foo', 'bar')])
        (args, kw) = mock_df.call_args
        self.assertListEqual([x.strip("\n").split("\t")
                              for x in self.data[8:]], args[0])
