# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import unittest
import logging
from snakemakelib.config import init_sml_config
from snakemakelib.bio.ngs.rnaseq.tuxedo import opt_read_group

logger = logging.getLogger(__name__)

def setUp():
    """Setup test_tuxedo test fixtures"""
    cfg = {
        'bio.ngs.settings' : {
            'run_id_re' : (("platform-unit", "date", "_", "sample"), "([0-9])_([0-9]+)_([A-Z0-9]+XX)_(P[0-9]+_[0-9]+)"), 
            'read_group_keys' : ("id", "sample", "library", "description", "platform-unit", "center", "date", "platform"),
            'center' : 'mycenter',
            'platform' : 'Illumina',
        },
    }
    init_sml_config(cfg)

def tearDown():
    """Teardown test_tuxedo test fixtures"""
    init_sml_config({})

class TestTuxedo(unittest.TestCase):
    """Test snakemakelib.bio.ngs.tuxedo module"""
    def test_opt_read_group(self):
        """Test generating read group"""
        s = opt_read_group("2_120924_AC003CCCXX_P001_102")
        self.assertEqual(s, "--rg-center mycenter --rg-date 2012-09-24 --rg-description 2_120924_AC003CCCXX_P001_102 --rg-id 2_120924_AC003CCCXX_P001_102 --rg-platform Illumina --rg-platform-unit 2 --rg-sample P001_102")
