# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
import os
import sys 
import unittest
import logging
from unittest.mock import patch
from nose.tools import raises
import pandas as pd
from snakemakelib.results import Results
from snakemakelib.exceptions import SamplesException, OutputFilesException

logging.basicConfig(level=logging.DEBUG)

class Foo(Results):
    _keys = ['foo', 'bar']
    def __init__(self, *args, **kw):
        super(Foo, self).__init__(*args, **kw)
        
    def _collect_results(self):
        pass

class TestResults(unittest.TestCase):
    """Test results base class"""
    def setUp(self):
        self.f = Foo(inputs=[("foo", "bar"), ("foofoo", "barbar")])

    @raises(NotImplementedError)
    def test_base_init(self):
        Results()

    @raises(SamplesException)
    def test_init_wrong_inputs(self):
        Results(inputs=["foo", "bar"])

    def test_init(self):
        Foo()

    def test_init_args(self):
        f = Foo(inputs=[("foo", "bar"), ("foofoo", "barbar")])

    @raises(OutputFilesException)
    def test_save_wrong_outputs(self):
        self.f.save([])

    @patch('pandas.DataFrame.to_csv')
    def test_save(self, mock_to_csv):
        self.f.save(["f", "b"])

    @patch('pandas.DataFrame')
    def test_parse1(self, mock_parse):
        self.f.parse_data([["foo", "bar"],["foofoo", "barbar"]])
        (args, _) = mock_parse.call_args
        self.assertEqual([["foo", "bar"],["foofoo", "barbar"]], args[0])

    @patch('pandas.DataFrame')
    def test_parse2(self, mock_parse):
        self.f.parse_data([["foo", "bar"],["foofoo", "barbar"]], rs=("foo", None))
        (args, _) = mock_parse.call_args
        self.assertEqual([["foo", "bar"], ["foofoo", "barbar"]], args[0])

    @patch('pandas.DataFrame')
    def test_parse3(self, mock_parse):
        self.f.parse_data([["foo", "bar"],["foofoo", "barbar"]], rs=(None, "barbar"))
        (args, _) = mock_parse.call_args
        self.assertEqual([["foo", "bar"]], args[0])


    @patch('pandas.DataFrame')
    def test_parse4(self, mock_parse):
        self.f.parse_data([["foo", "bar"],["foofoo", "barbar"]], rs=(None, None), skip=1)
        (args, _) = mock_parse.call_args
        self.assertEqual([["foofoo", "barbar"]], args[0])


    @patch('pandas.read_csv')
    def test_load_data_frame(self, mock_read_csv):
        f = Foo(inputs=[("foo", "bar"), ("foofoo", "barbar")])
        f.load_data_frame("foo")
        (args, kw) = mock_read_csv.call_args
        self.assertEqual(('foo',), args)
