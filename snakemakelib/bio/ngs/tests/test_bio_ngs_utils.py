# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import os
from unittest.mock import patch
from snakemakelib.bio.ngs.utils import find_files

class TestFindFiles(unittest.TestCase):
    def setUp(self):
        """Setup text fixtures"""
        self.walk = [
            [os.curdir, ['foo', 'bar'], ['1_121023_FLOWCELL_FOO.fastq.gz', 'bar.txt']],
            ['./foo', [], ['foo.txt', '1_121023_FLOWCELL_BAR.fastq.gz']],
            ['./bar', [], ['bar.txt']],
            ]

    @patch('snakemakelib.bio.ngs.utils.os.walk')
    def test_find_fastq_files(self, mock_walk):
        """Find fastq files using match"""
        mock_walk.return_value = self.walk
        f = find_files(regexp="\w+.fastq.gz")
        self.assertListEqual(f, ['./1_121023_FLOWCELL_FOO.fastq.gz', './foo/1_121023_FLOWCELL_BAR.fastq.gz'])

    @patch('snakemakelib.bio.ngs.utils.os.walk')
    def test_find_files_search(self, mock_walk):
        """Find files using search"""
        mock_walk.return_value = self.walk
        f = find_files(regexp=".fastq", search=True)
        self.assertListEqual(f, ['./1_121023_FLOWCELL_FOO.fastq.gz', './foo/1_121023_FLOWCELL_BAR.fastq.gz'])
