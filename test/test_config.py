# Copyright (C) 2014 by Per Unneberg
import os
import sys 
import unittest
import logging
from nose.tools import raises
from snakemakelib.config import sml_config, BaseConfig, update_sml_config
from snakemake.logging import logger

logging.basicConfig(level=logging.DEBUG)

def fun1():
    return "fun1"

class TestBaseConfig(unittest.TestCase):
    def setUp(self):
        self.cfg = BaseConfig({'foo':'bar'})

    def tearDown(self):
        del self.cfg

    def test_create_cfg_from_dict(self):
        """Test create configuration from dictionary"""
        self.assertIsInstance(self.cfg, BaseConfig)
        self.assertDictEqual(self.cfg, {'foo':'bar'})
        self.assertListEqual(self.cfg.sections, ['foo'])

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
        del cfg

    def test_create_cfg_from_args_kwargs(self):
        """Test create configuration from *args and **kwargs"""
        cfg = BaseConfig(foo="bar", bar=BaseConfig(foo="bar", bar="foo"), **{'foobar':'bar', 'barfoo':BaseConfig({'foo':'bar'})})
        self.assertSetEqual (set(cfg.sections), set(['bar', 'foo', 'foobar', 'barfoo']))
        self.assertIsInstance(cfg['bar'], BaseConfig)
        self.assertIsInstance(cfg['barfoo'], BaseConfig)
        del cfg

    @raises(TypeError)
    def test_create_cfg_from_nested_dict_with_subsection_dict(self):
        """Test create configuration from dictionary in which a subsection is a dict. Should raise TypeError error."""
        cfg = BaseConfig({'foo':BaseConfig({'bar':{'foo':'bar'}})})
        del cfg

    @raises(TypeError)
    def test_add_section_dict(self):
        """Test adding a section to config as dict"""
        self.cfg.add_section({'foobar':'bar'})

    def test_add_section_str(self):
        """Test adding a section to config as str"""
        self.cfg.add_section('foobar')
        self.assertSetEqual(set(self.cfg.sections), set(['foo', 'foobar']))


    def test_update_config(self):
        """Test updating configuration with another configuration object. """
        cfg2 = BaseConfig({'bar':'foo'})
        self.cfg.update(cfg2)
        self.assertDictEqual(self.cfg, {'foo':'bar', 'bar':'foo'})
        del cfg2

    def test_update_config_same_section(self):
        """Test updating configuration with another configuration objet whose sections overlap. Note that this will overwrite the first configuration. FIXME: should this be intended behaviour?"""
        cfg2 = BaseConfig({'foo':'foobar'})
        self.cfg.update(cfg2)
        self.assertDictEqual(self.cfg, {'foo':'foobar'})
        del cfg2

    def test_update_config_nested(self):
        """Test updating configuration with another nested configuration"""
        cfg2 = BaseConfig({'bar': BaseConfig({'foo':'bar'})})
        self.cfg.update(cfg2)
        self.assertDictEqual(self.cfg, {'foo':'bar', 'bar':{'foo':'bar'}})
        del cfg2


    @raises(TypeError)
    def test_update_config_nested_dict(self):
        """Test updating configuration with another nested configuration that has a dict section"""
        cfg2 = BaseConfig({'bar': {'foo':'bar'}})
        self.cfg.update(cfg2)
        del cfg2

    @raises(TypeError)
    def test_setting_config_section_to_dict(self):
        """Test setting a configuration section to a dict"""
        self.cfg['foo'] = {'foo':'bar'}

class TestSmlConfig(unittest.TestCase):
    def setUp(self):
        self.cfg = BaseConfig({
            'bar' : BaseConfig({
                'foo' : 'customfoo',
            })
        })
        self.cfg_nested = BaseConfig({
            'bar' : BaseConfig({
                'bar' : BaseConfig({
                    'bar':'customfoo'
                })
            })
        })

        self.default = BaseConfig({
            'foo':'bar',
            'bar' : BaseConfig({
                'foo' : 'foobar',
                'bar' : 'foo'
            })
        })
        self.default_nested = BaseConfig({
            'foo':'bar',
            'bar' : BaseConfig({
                'foo' : 'foobar',
                'bar' : BaseConfig({
                    # bar key is missing so should not be possible to do both nested?!?
                    'foo' : 'bar'
                })
            })
        })
        
    def tearDown(self):
        del self.cfg
        del self.default

    # def test_sml_config(self):
    #     """Test printing global sml config"""
    #     print ("Global configuration " , sml_config)


    def test_update_sml_config_with_default(self):
        """Test updating a configuration object skipping values that are already in use. Overriding dict.update will not do as its original intedend behaviour is needed.

        What is the desired behaviour? 

        1. In Snakefile user modifies a custom configuration object
        2. Relevant include files are loaded with default settings.
        3. Default settings need to be updated with custom config at once so that custom changes are reflected in rules (is this true?)
        """
        self.cfg = update_sml_config(self.cfg, self.default)
        self.assertDictEqual(self.cfg, {'bar': {'bar': 'foo', 'foo': 'customfoo'}, 'foo': 'bar'})


    def test_update_sml_config_with_default_nested(self):
        self.cfg = update_sml_config(self.cfg, self.default_nested)
        self.assertDictEqual(self.cfg, {'foo': 'bar', 'bar': {'foo': 'customfoo', 'bar': {'foo': 'bar'}}})


    def test_update_sml_config_with_both_nested(self):
        """Test updating a configuration object where both are nested. Note that in this example  self.cfg_nested has a key (section) not defined in default so a warning should be raised. In other words, at a given level, if default is a BaseConfig, the keys in config should be a subset of keys in default."""
        self.cfg_nested = update_sml_config(self.cfg_nested, self.default_nested)
        self.assertDictEqual(self.cfg_nested, {'foo': 'bar', 'bar': {'foo': 'foobar', 'bar': {'foo': 'bar', 'bar': 'customfoo'}}})

    @raises(TypeError)
    def test_update_sml_config_with_cfg_nested_missing_bc(self):
        cfg2 = BaseConfig({
            'bar' : BaseConfig({
                'bar' : 'test'
                })
            })
        cfg2 = update_sml_config(cfg2, self.default_nested)
