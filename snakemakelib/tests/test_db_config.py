# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
import os
import unittest
import logging
from nose.tools import raises
from snakemakelib.bio.ngs.db import ref, index

logging.basicConfig(level=logging.DEBUG)

class TestDbConfig(unittest.TestCase):
    """Test db configuration utilities"""
    def setUp(self):
        def ref():
            return os.path.abspath(os.path.join(os.curdir, "seq", "reference.fa"))
        self.ref = ref

    def test_ref_wo_build(self):
        """Test reference function return value when build is missing"""
        cfg = {'bio.ngs.settings' : {'db' : {'build' : None, 'ref':"foo/bar"}}}
        self.assertEqual(ref("foobar", cfg['bio.ngs.settings']['db']), "foo/foobar")

    def test_ref_wo_build_config(self):
        """Test reference function return value when build present but build_config is not"""
        cfg = {'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None, 'ref' : "foo/bar"}}}
        self.assertEqual(ref("foobar", cfg['bio.ngs.settings']['db']), "foo/foobar")

    def test_custom_ref(self):
        """Test setting a custom ref"""
        config = {'bio.ngs.settings' : {'db' : {'ref' : self.ref()}}}
        cfg = config['bio.ngs.settings']
        cfg_d = dict(cfg['db'])
        self.assertIsInstance(cfg['db']['ref'], str)
        self.assertEqual(cfg['db']['ref'],
                         os.path.abspath(
                             os.path.join(os.curdir, "seq", "reference.fa")))

    @raises(TypeError)
    def test_index_missing_application(self):
        """Test calling index without applying application"""
        config = {'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None, 'ref' : self.ref}}, 'bio.ngs.align.bwa': {'index' : index}}
        bwa_cfg = config['bio.ngs.align.bwa']
        index(ref= config['bio.ngs.settings']['db']['ref'], index=bwa_cfg['index'])
