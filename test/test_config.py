# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import unittest
import logging
from nose.tools import raises
from snakemakelib.config import sml_config, BaseConfig
from snakemake.logging import logger

logging.basicConfig(level=logging.DEBUG)

def fun1():
    return "fun1"



class TestBaseConfig(unittest.TestCase):
    def setUp(self):
        self.cfg_foo = BaseConfig({'foo':'bar'})
        self.cfg_bar = BaseConfig({'bar':'foo', 'fun':fun1()})

    def test_create_cfg_from_dict(self):
        """Test create configuration from dictionary"""
        cfg = BaseConfig({'foo':'bar'})
        self.assertIsInstance(cfg, BaseConfig)
        self.assertDictEqual(cfg, {'foo':'bar'})
        self.assertListEqual(cfg.sections, ['foo'])

    def test_create_cfg_from_nested_dict(self):
        """Test create configuration from nested dictionary"""
        cfg = BaseConfig({'foo': BaseConfig({'bar':'foobar'})})
        self.assertIsInstance(cfg, BaseConfig)
        self.assertIsInstance(cfg['foo'], BaseConfig)
        self.assertDictEqual(cfg['foo'], {'bar':'foobar'})
        self.assertListEqual(cfg.sections, ['foo'])
        self.assertListEqual(cfg['foo'].sections, ['bar'])

    def test_create_cfg_from_args(self):
        """Test create configuration from *args"""
        cfg = BaseConfig(foo="bar", bar=BaseConfig(foo="bar", bar="foo"))
        self.assertSetEqual (set(cfg.sections), set(['bar', 'foo']))
        self.assertIsInstance(cfg['bar'], BaseConfig)

    def test_create_cfg_from_args_kwargs(self):
        """Test create configuration from *args and **kwargs"""
        cfg = BaseConfig(foo="bar", bar=BaseConfig(foo="bar", bar="foo"), **{'foobar':'bar', 'barfoo':BaseConfig({'foo':'bar'})})
        self.assertSetEqual (set(cfg.sections), set(['bar', 'foo', 'foobar', 'barfoo']))
        self.assertIsInstance(cfg['bar'], BaseConfig)
        self.assertIsInstance(cfg['barfoo'], BaseConfig)

    @raises(TypeError)
    def test_create_cfg_from_nested_dict_with_subsection_dict(self):
        """Test create configuration from dictionary in which a subsection is a dict. Should raise TypeError error."""
        cfg = BaseConfig({'foo':BaseConfig({'bar':{'foo':'bar'}})})

    @raises(TypeError)
    def test_add_section_dict(self):
        """Test adding a section to config as dict"""
        self.cfg_foo.add_section({'foobar':'bar'})

    def test_add_section_str(self):
        """Test adding a section to config as str"""
        self.cfg_foo.add_section('foobar')

    def test_update_config(self):
        """Test updating configuration with another configuration objet."""
        # The question is whether to override dict.update or have a
        # separate function?
        pass


    @raises(TypeError)
    def test_setting_config_section_to_dict(self):
        """Test setting a configuration section to a dict"""
        cfg = BaseConfig({'foo':'bar'})
        cfg['foo'] = {'foo':'bar'}

class TestGlobalConfig(unittest.TestCase):
    def test_sml_config(self):
        """Test printing global sml config"""
        print ("Global configuration " , sml_config)

