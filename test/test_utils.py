# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
from snakemakelib.utils import isoformat

class TestUtils(unittest.TestCase):
    """Test snakemakelib.utils functions"""
    def test_isoformat(self):
        """Test isoformatting function"""
        s = "120924"
        self.assertEqual(isoformat(s), "2012-09-24")

