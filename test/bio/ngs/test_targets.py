# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import os
import unittest
import logging
import csv
import io
from unittest.mock import patch
from nose.tools import raises
import snakemakelib.bio.ngs.targets
from snakemakelib.bio.ngs.targets import generic_target_generator, generic_target_generator2
from snakemakelib.bio.ngs.regexp import ReadGroup

logger = logging.getLogger(__name__)

class TestTargetGenerator(unittest.TestCase):
    @raises(AssertionError)
    def test_regexp_call(self):
        """Make sure called regexp is of class RegexpDict"""
        generic_target_generator(tgt_re="regexp", **{})

    def test_from_command_line(self):
        """Generate targets from command line options. Match foo/1_BAR123XX_foo"""
        tgts = generic_target_generator(tgt_re=ReadGroup(os.path.join("(?P<SM>[a-z]+)", "(?P<PU1>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P=SM)")), **{'samples' : ['foo', 'bar'], 'runs' : ['1_BAR123XX_foo', '2_BAR456XX_bar'], 'target_suffix' : ".bam"})
        self.assertListEqual(tgts, ['bar/2_BAR456XX_bar.bam', 'foo/1_BAR123XX_foo.bam'])

    @patch('snakemakelib.bio.ngs.targets.smllogger.warn')
    def test_from_command_line_unequal(self, mock_warn):
        """Generate targets from command line options; len(samples) != len(runs)"""
        generic_target_generator(tgt_re=ReadGroup(""), **{'samples':['foo', 'bar'], 'runs':['foo']})
        self.assertTrue(mock_warn.called)

    @patch('snakemakelib.bio.ngs.targets.os.path.exists')
    @patch('builtins.open')
    def test_from_samplesheet(self, mock_open, mock_exists):
        """Generate targets from samplesheet information"""
        mock_exists.return_value = True
        mock_open.return_value = io.StringIO("\n".join(["FCID,Lane,SampleID", "FC1,1,S1", "FC2,2,S2", "FC3,1,S3"]))
        tgts = generic_target_generator(tgt_re=ReadGroup("(?P<PU1>[0-9]+)_(?P<PU2>[A-Z0-9]+)_(?P<SM>[0-9A-Z]+)"), **{'sampleinfo':'foo', 'sample_column_map': {'FCID':'PU2', 'Lane':'PU1', 'SampleID':'SM'}})
        self.assertListEqual(tgts, ['1_FC1_S1', '1_FC3_S3', '2_FC2_S2'])
        
    @patch('snakemakelib.bio.ngs.targets.os.path.exists')
    @patch('builtins.open')
    def test_from_samplesheet_subsample(self, mock_open, mock_exists):
        """Generate targets from samplesheet information, selecting one sample via 'samples'"""
        mock_exists.return_value = True
        mock_open.return_value = io.StringIO("\n".join(["FCID,Lane,SampleID", "FC1,1,S1", "FC2,2,S2", "FC3,1,S3"]))
        tgts = generic_target_generator(tgt_re=ReadGroup("(?P<PU1>[0-9]+)_(?P<PU2>[A-Z0-9]+)_(?P<SM>[0-9A-Z]+)"), **{'sampleinfo':'foo', 'sample_column_map': {'FCID':'PU2', 'Lane':'PU1', 'SampleID':'SM'}, 'samples':['S1']})
        self.assertListEqual(tgts, ['1_FC1_S1'])

    @patch('snakemakelib.bio.ngs.targets.find_files')
    def test_from_input_files(self, mock_ff):
        """Generate targets from input file information"""
        mock_ff.return_value = ['S1/1_FC1_S1_1.fastq.gz', 'S2/2_FC2_S2_1.fastq.gz']
        tgts = generic_target_generator(tgt_re=ReadGroup("(?P<SM>[0-9A-Z]+)/(?P<PU1>[0-9]+)_(?P<PU2>[A-Z0-9]+)_(?P=SM)"), **{'target_suffix':'.sort.merge.bam'})
        self.assertListEqual(tgts, ['S1/1_FC1_S1.sort.merge.bam', 'S2/2_FC2_S2.sort.merge.bam'])

    @patch('snakemakelib.bio.ngs.targets.find_files')
    def test_from_input_files_src_re(self, mock_ff):
        """Generate targets from input file information; use different src_re"""
        mock_ff.return_value = ['S1/1_FC1_S1_1.fastq.gz', 'S2/2_FC2_S2_1.fastq.gz']
        tgts = generic_target_generator(tgt_re=ReadGroup("(?P<SM>[0-9A-Z]+)/(?P=SM)"), src_re=ReadGroup("(?P<SM>[0-9A-Z]+)/(?P<PU1>[0-9]+)_(?P<PU2>[A-Z0-9]+)_(?P=SM)"), **{'target_suffix':'.sort.merge.bam'})
        self.assertListEqual(tgts, ['S1/S1.sort.merge.bam', 'S2/S2.sort.merge.bam'])

