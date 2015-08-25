# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
import unittest
import logging
import io
import mock
from nose.tools import raises
from snakemakelib.bio.ngs.tools.sratools import register_metadata

logger = logging.getLogger(__name__)

class TestSraTools(unittest.TestCase):
    """Test sratools functionality"""
    def setUp(self):
        self.metadata = io.StringIO("\n".join(["SampleName,Run",
                                               "Sample1,Run1S1",
                                               "Sample1,Run2S1",
                                               "Sample2,Run1S2"]))

    @raises(Exception)
    def test_register_metadata(self):
        """Test registering metadata for non-existent file"""
        register_metadata("foo.csv", config = {})

    @mock.patch('snakemakelib.bio.ngs.tools.sratools.update_config')
    @mock.patch('builtins.open')
    def test_read(self, mock_open, mock_update):
        """Test reading mock indata and check update arguments"""
        mock_open.return_value = self.metadata
        cfg = register_metadata("foo/bar", config = {})
        args, kw = mock_update.call_args
        self.assertDictEqual(args[0]['bio.ngs.tools.sratools']['_run2sample'], {'Run2S1': 'Sample1', 'Run1S2': 'Sample2', 'Run1S1': 'Sample1'})
        self.assertEqual(args[0]['bio.ngs.tools.sratools']['_datadir'], 'foo')
        
