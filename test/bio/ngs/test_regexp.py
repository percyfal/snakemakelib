# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import re
import os
from unittest.mock import patch
from itertools import groupby
from nose.tools import raises
from snakemakelib.bio.ngs.regexp import RegexpDict, SampleRegexp, ReadGroup, DisallowedKeyException, MissingRequiredKeyException

class TestRegexpDict(unittest.TestCase):
    @raises(TypeError)
    def test_regexpdict_init_empty(self):
        rd = RegexpDict()

    @raises(DisallowedKeyException)
    def test_regexpdict_init_wrong_keys(self):
        rd = RegexpDict("(?P<PU1>[A-Z0-9]+)_(?P<PU2>[0-9]+)")

    def test_regexpdict_basename_pattern(self):
        rd = SampleRegexp(os.path.join("(?P<SM>[a-z]+)", "(?P<PU>[0-9]+)_(?P=SM)"))
        self.assertEqual(rd.pattern, "(?P<SM>[a-z]+)/(?P<PU>[0-9]+)_(?P=SM)")
        self.assertEqual(rd.basename_pattern, "(?P<PU>[0-9]+)_(?P<SM>[a-z]+)")

    def test_regexpdict_relative_basename_pattern(self):
        """Test basename of sample regexp without leading paths"""
        rd = SampleRegexp(os.path.join("(?P<PU>[0-9]+)_(?P<SM>[a-z]+)"))
        self.assertEqual(rd.pattern, "(?P<PU>[0-9]+)_(?P<SM>[a-z]+)")
        self.assertEqual(rd.basename_pattern, "(?P<PU>[0-9]+)_(?P<SM>[a-z]+)")
        
class TestSampleRegexp(unittest.TestCase):
    @raises(TypeError)
    def test_sampleregexp_init_empty(self):
        rd = SampleRegexp()

    @raises(MissingRequiredKeyException)
    def test_sampleregexp_missing_required_key(self):
        rd = SampleRegexp("(?:[A-Z0-9]+)_(?P<PATH>[0-9]+)")
    
    def test_sampleregexp_parse_unable(self):
        rd = SampleRegexp("(?:[A-Z0-9]+)_(?P<SM>[A-Za-z0-9]+)")
        self.assertDictEqual(rd.parse("bar/foo_007"), {})

    @raises(DisallowedKeyException)
    def test_sampleregexp_wrong_indexed_key(self):
        rd = SampleRegexp("(?:[A-Z0-9]+)_(?P<SM>[0-9]+)_(?P<ID1>[0-9]+)")

    def test_sampleregexp_correct_indexed_key(self):
        rd = SampleRegexp(os.path.join("(?P<SM>[A-Za-z0-9]+)", "(?:[A-Z0-9]+)_(?P<PU2>[0-9]+)_(?P<PU1>[0-9]+)"))
        self.assertEqual(rd.fmt, "{SM}/{PU2}_{PU1}")
        rd.parse("BAR/FOO_007_1")
        self.assertDictEqual(rd, {'PATH': 'BAR', 'PU': '1_007', 'PU2': '007', 'PU1': '1', 'SM' : 'BAR'})

class TestGroupBy(unittest.TestCase):
    def test_list(self):
        r = re.compile("(?P<SM>[A-Z0-9]+)_(?P<ID2>[0-9]+)_(?P<ID1>[0-9]+)")
        keymap = sorted([(re.sub("[0-9]+$", "", k), k)  if re.search("[0-9]+$", k) else (k, k) for k in list(r.groupindex.keys())])
        _keymap = {k:[y[1] for y in list(v)] for (k,v) in groupby(keymap, lambda x: x[0])}
        self.assertDictEqual(_keymap, {'ID': ['ID1', 'ID2'], 'SM': ['SM']})

    def test_make_format(self):
        """Covert named group to format tags.

        From (?P<SM>[A-Z0-9]+)/(?P=SM)_(?P<ID2>[0-9]+)_(?P<ID1>[0-9]+) generate {SM}/{SM}_{ID2}_{ID1}"""
        r = re.compile(os.path.join("(?P<SM>[A-Z0-9]+)", "(?P=SM)_(?P<ID2>[0-9]+)_(?P<ID1>[0-9]+)"))
        m = re.findall("(\(\?P[<=](\w+)>?|({sep}))".format(sep=os.sep), r.pattern)
        fmt = re.sub("_{sep}_".format(sep=os.sep), os.sep, ("_".join("{" + x[1] + "}"  if x[1] else x[2] for x in m)))
        self.assertEqual(fmt, "{SM}/{SM}_{ID2}_{ID1}")

class TestReadGroup(unittest.TestCase):
    """Test ReadGroup class"""
    def test_rg_init(self):
        """Test initializing ReadGroup"""
        rg = ReadGroup("(?P<ID>[A-Za-z0-9]+)", ID='test', DT="120924")
        self.assertEqual(str(rg), '--date 2012-09-24 --identifier test')
        rg = ReadGroup("(?P<ID>[A-Za-z0-9]+)", **{'ID':'test', 'DT':"120924"})
        self.assertEqual(str(rg), '--date 2012-09-24 --identifier test')
        
    def test_rg_parse_illumina_like(self):
        """Test parsing illumina-like-based file names"""
        rg = ReadGroup("(?P<SM>[a-zA-Z0-9_]+)/(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)/(?P<PU1>[0-9])_(?P=DT)_(?P=PU2)_(?P=SM)")
        rg.parse("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz")
        self.assertEqual("--date 2012-10-15 --identifier 1_121015_BB002BBBXX_P001_101 --platform-unit 1_BB002BBBXX --sample P001_101", str(rg))

    def test_rg_fn(self):
        """Test initializing read group class and setting function"""
        rg_fn = ReadGroup("(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)").parse
        s = rg_fn("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz")
        self.assertEqual("--date 2012-10-15 --identifier 1_121015_BB002BBBXX_P001_101 --platform-unit 1_BB002BBBXX --sample P001_101", str(s))

    @raises(DisallowedKeyException)
    def test_rg_disallowed_key(self):
        """Test setting a read group object with a key not present in allowed keys"""
        rg = ReadGroup("(?P<PU1>[0-9])_(?P<DATE>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P<SM>P[0-9]+_[0-9]+)")
        s = (rg.parse("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz"))

# NB: all regexp names must be relative, and complete, to the working
# directory. IOW, the target generator must return paths of this
# format. Still, test full path names below just to make sure we catch
# anomalous cases.
class TestParseFunctionality(unittest.TestCase):
    def setUp(self):
        self.re = r"(?P<SM>\w+)/(?P<PU>[A-Za-z0-9_]+)/(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P=SM)"
        self.full_re = r"(?:[\.\w\/]+)?\/(?P<SM>\w+)/(?P<PU>[A-Za-z0-9_]+)/(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P=SM)"
        self.fn = "P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz"
        self.full_fn = "../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz"

    def tearDown(self):
        del self.re, self.full_re, self.fn, self.full_fn

    def test_relpath(self):
        m = re.search(self.re, self.fn)
        self.assertDictEqual(m.groupdict(), {'SM': 'P001_101', 'DT': '121015', 'PU1': '1', 'PU2': 'BB002BBBXX', 'PU': '121015_BB002BBBXX'})

    def test_curdir_path(self):
        m = re.search(self.full_re, os.path.join(os.curdir, self.fn))
        self.assertDictEqual(m.groupdict(), {'SM': 'P001_101', 'DT': '121015', 'PU1': '1', 'PU2': 'BB002BBBXX', 'PU': '121015_BB002BBBXX'}) 
        
    def test_pardir_path(self):
        m = re.search(self.full_re, self.full_fn)
        self.assertDictEqual(m.groupdict(), {'SM': 'P001_101', 'DT': '121015', 'PU1': '1', 'PU2': 'BB002BBBXX', 'PU': '121015_BB002BBBXX'})
        
    def test_abspath(self):
        m = re.search(self.full_re, os.path.abspath(self.full_fn))
        self.assertDictEqual(m.groupdict(), {'SM': 'P001_101', 'DT': '121015', 'PU1': '1', 'PU2': 'BB002BBBXX', 'PU': '121015_BB002BBBXX'})

    @patch('snakemakelib.bio.ngs.regexp.re.match')
    def test_pardir(self, mock_re):
        mock_re.return_value = None
        rg = ReadGroup("(?P<SM>[a-zA-Z0-9]+)/(?P<PU>[A-Za-z0-9]+)/(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P=SM)")
        rg.parse("../data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz", "")
        (args, kw) = mock_re.call_args
        self.assertTrue(args[0].startswith('(?:[\\.\\w\\/]+)?\\/'))

    @patch('snakemakelib.bio.ngs.regexp.re.match')
    def test_curdir(self, mock_re):
        mock_re.return_value = None
        rg = ReadGroup(r"(?P<SM>[a-zA-Z0-9]+)/(?P<PU>[A-Za-z0-9]+)/(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P=SM)")
        rg.parse("./data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz", "")
        (args, kw) = mock_re.call_args
        self.assertTrue(args[0].startswith('(?:[\\.\\w\\/]+)?\\/'))

    @patch('snakemakelib.bio.ngs.regexp.re.match')
    def test_os_sep(self, mock_re):
        mock_re.return_value = None
        rg = ReadGroup("(?P<SM>[a-zA-Z0-9]+)/(?P<PU>[A-Za-z0-9]+)/(?P<PU1>[0-9])_(?P<DT>[0-9]+)_(?P<PU2>[A-Z0-9]+XX)_(?P=SM)")
        rg.parse("/data/projects/J.Doe_00_01/P001_101/121015_BB002BBBXX/1_121015_BB002BBBXX_P001_101_1.fastq.gz", "")
        (args, kw) = mock_re.call_args
        self.assertTrue(args[0].startswith('(?:[\\.\\w\\/]+)?\\/'))
