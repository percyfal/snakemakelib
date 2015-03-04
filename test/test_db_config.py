# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
import os
import sys 
import unittest
import logging
from nose.tools import raises
from snakemakelib.config import init_sml_config, update_sml_config, get_sml_config
from snakemakelib.bio.ngs.db import ref, index

logging.basicConfig(level=logging.DEBUG)

class TestDbConfig(unittest.TestCase):
    """Test db configuration utilities"""
    def setUp(self):
        def ref():
            return os.path.abspath(os.path.join(os.curdir, "seq", "reference.fa"))
        self.ref = ref

    @raises(ValueError)
    def test_ref_wo_build(self):
        """Test reference function return value when build is missing"""
        init_sml_config({'bio.ngs.settings' : {'db' : {'build' : None}}})
        ref()

    @raises(ValueError)
    def test_ref_wo_build_config(self):
        """Test reference function return value when build present but build_config is not"""
        init_sml_config({'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None}}})
        ref()

    def test_custom_ref(self):
        """Test setting a custom ref"""
        init_sml_config({'bio.ngs.settings' : {'db' : {'ref' : self.ref()}}})
        cfg = get_sml_config('bio.ngs.settings')
        cfg_d = dict(cfg['db'])
        self.assertIsInstance(cfg['db']['ref'], str)
        self.assertEqual(cfg['db']['ref'],
                         os.path.abspath(
                             os.path.join(os.curdir, "seq", "reference.fa")))

    def test_custom_ref_fn(self):
        """Test setting a custom ref function"""
        init_sml_config({'bio.ngs.settings' : {'db' : {'ref' : self.ref}}})
        cfg = get_sml_config('bio.ngs.settings')
        cfg_d = dict(cfg['db'])
        self.assertEqual(str(type(cfg_d['ref'])), "<class 'function'>")
        self.assertIsInstance(cfg['db']['ref'], str)
        self.assertEqual(cfg['db']['ref'], os.path.abspath(
            os.path.join(os.curdir, "seq", "reference.fa")))

    def test_index(self):
        """Test getting the index for a given application"""
        init_sml_config({'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None, 'ref' : self.ref}}, 'bio.ngs.align.bwa': {'index' : index}})
        cfg = get_sml_config('bio.ngs.settings')
        bwa_cfg = get_sml_config('bio.ngs.align.bwa')
        self.assertEqual(bwa_cfg['index', 'bwa'], os.path.abspath(
            os.path.join(os.curdir, "bwa", "reference")))

    @raises(TypeError)
    def test_index_missing_application(self):
        """Test calling index without applying application"""
        init_sml_config({'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None, 'ref' : self.ref}}, 'bio.ngs.align.bwa': {'index' : index}})
        bwa_cfg = get_sml_config('bio.ngs.align.bwa')
        bwa_cfg['index']
