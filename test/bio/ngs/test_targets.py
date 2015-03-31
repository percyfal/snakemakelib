# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import os
import unittest
import logging
import csv
from nose.tools import raises
from snakemakelib.bio.ngs.targets import generic_target_generator
from snakemakelib.bio.ngs.regexp import ReadGroup

logger = logging.getLogger(__name__)

class TestTargetGenerator(unittest.TestCase):
    """Test target generating functions and utilities"""
    def setUp(self):
        self._sampleinfo = ["FCID,Lane,SampleID", "FC1,1,S1", "FC2,2,S2", "FC3,1,S3"]
        def _parse_sampleinfo(fn):
            reader = csv.DictReader(fn)
            keymap = {'FCID':'PU2', 'Lane':'PU1', 'SampleID':'SM'}
            reader.fieldnames = [keymap[f] for f in reader.fieldnames]
            return reader
        
        self.generic_cfg = {
            'run_id_re' : "(?P<PU1>[0-9])_(?P<PU2>[A-Z0-9]+)_(?P<SM>[A-Z0-9]+)",
            'run_id_pfx_re' : "(?P<PU1>[0-9])_(?P<PU2>[A-Z0-9]+)_(?P<SM>[A-Z0-9]+)",
            'run_id_pfx_fmt' : os.path.join("{SM}", "{PU2}", "{PU1}_{PU2}_{SM}"), 
            'sample_pfx_fmt' : os.path.join("{SM}", "{SM}"),
            'samples':["S1", "S2", "S3"],
            'runs' : ["1_FC1_S1", "2_FC2_S2", "1_FC3_S3"],
            'sampleinfo' : '',
            }
        self.cfg = {
            'run_id_re' : "(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)",
            'run_id_pfx_re' : "(?P<PU1>[0-9])_(?P<PU2>[0-9]+_[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)",
            'run_id_pfx_fmt' : os.path.join("{SM}", "{PU2}", "{PU1}_{PU2}_{SM}"), 
            'sample_pfx_fmt' : os.path.join("{SM}", "{SM}"),
            'platform_unit_fn' : lambda x: {'PU1':x[2], 'PU2':x[1]},
            'sample_column_name' : 'SM',
            'run_column_name' : 'PU',
        }
        self.cfg2 = {
            'samples':["S2"],
            'runs' : [],
            'sampleinfo' : '',
        }
        self.cfg2.update(self.cfg)
        self.cfg3 = {
            'samples' : ['S2'],
            'runs' : [],
            'sampleinfo' : _parse_sampleinfo(self._sampleinfo),
            }
        self.cfg3.update(self.cfg)
        self.cfg4 = {
            'samples' : [],
            'runs' : [],
            'sampleinfo' : _parse_sampleinfo(self._sampleinfo),
            }
        self.cfg4.update(self.cfg)
        self.cfg5 = {
            'samples' : [],
            'runs' : [],
            'sampleinfo' : '',
            }
        self.cfg5.update(self.cfg)
        self.path = "../data/projects/J.Doe_00_01"

    # Order of preference:
    # 1. cfg values
    # 2. samples
    # 3. generate name from input files (generated from run_id_re)
    # 4. possibly search for SampleSheet.csv, parse and generate names
    def test_from_generic_input(self):
        """Test generating targets from input files"""
        l = generic_target_generator(fmt=self.cfg['run_id_pfx_fmt'] + ".bam", rg=ReadGroup(regexp=self.cfg["run_id_pfx_re"]), cfg=self.cfg, path=self.path)
        self.assertEqual(l[0], '../data/projects/J.Doe_00_01/S1/FC1/1_FC1_S1.bam')

    @raises(Exception)
    def test_sample_only(self):
        """Test generating targets from sample only"""
        l = generic_target_generator(fmt=self.cfg['run_id_pfx_fmt'] + ".bam", rg=ReadGroup(regexp=self.cfg["run_id_pfx_re"]), cfg=self.cfg2, path=self.path)

    def test_sample_w_sampleinfo(self):
        """Test generating targets from sample only, with sampleinfo present.
        """
        l = generic_target_generator(fmt=self.cfg['run_id_pfx_fmt'] + ".bam", rg=ReadGroup(regexp=self.cfg["run_id_pfx_re"]), cfg=self.cfg3, path=self.path)
        self.assertEqual(len(l), 1)
        self.assertEqual(l[0], '../data/projects/J.Doe_00_01/S2/FC2/2_FC2_S2.bam')

    def test_sampleinfo(self):
        """Test generating targets from sampleinfo.
        """
        l = generic_target_generator(fmt=self.cfg['run_id_pfx_fmt'] + ".bam", rg=ReadGroup(regexp=self.cfg["run_id_pfx_re"]), cfg=self.cfg4, path=self.path)
        self.assertEqual(l[1], '../data/projects/J.Doe_00_01/S2/FC2/2_FC2_S2.bam')
        
    def test_from_input(self):
        """Test generating targets from input file names"""
        l = generic_target_generator(fmt=self.cfg['run_id_pfx_fmt'] + ".bam", rg=ReadGroup(regexp=self.cfg["run_id_pfx_re"] + "_1.fastq.gz$"), cfg=self.cfg5, path=self.path)
        self.assertEqual(l[1], '../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101.bam')
