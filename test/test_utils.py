# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
from snakemakelib.stat import is_installed
from snakemakelib.utils import isoformat

class TestUtils(unittest.TestCase):
    """Test snakemakelib.utils functions"""
    def test_isoformat(self):
        """Test isoformatting function"""
        s = "120924"
        self.assertEqual(isoformat(s), "2012-09-24")

class TestStat(unittest.TestCase):
    """Test snakemakelib.stat functions"""
    def test_is_installed(self):
        """Test function for checking that program is installed or that path exists"""
        self.assertTrue (is_installed("ls"))
        self.assertTrue(is_installed("/bin/ls"))
        self.assertTrue(is_installed("/bin"))
        self.assertFalse(is_installed("ls2"))
