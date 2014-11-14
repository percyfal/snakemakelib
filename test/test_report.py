# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import unittest
import logging
from nose.tools import raises
from snakemake.report.picard import 
from snakemake.logging import logger

# FIXME: check in metrics files
class TestMergeDict(unittest.TestCase):
    def test_merge_two_dicts(self):
        """Test merging two dictionaries where we know we have the same columns"""
        
