# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import logging
import re
import os
from unittest.mock import patch
from nose.tools import raises
from snakemakelib.bio.ngs.qc.qualimap import Qualimap

class TestQualimap(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        Qualimap()
