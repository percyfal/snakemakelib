# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import unittest
import logging
from snakemakelib.bio.ngs.rnaseq.tuxedo import TuxedoReadGroup

logger = logging.getLogger(__name__)

class TestTuxedoReadGroup(unittest.TestCase):
    """Test TuxedoReadGroup class"""
    def test_rg_init(self):
        """Test initializing TuxedoReadGroup"""
        rg = TuxedoReadGroup(regexp="test", ID='test', DT="120924")
        self.assertEqual(str(rg), '--rg-date 2012-09-24 --rg-id test')
        rg = TuxedoReadGroup(regexp="test", **{'ID':'test', 'DT':"120924"})
        self.assertEqual(str(rg), '--rg-date 2012-09-24 --rg-id test')

