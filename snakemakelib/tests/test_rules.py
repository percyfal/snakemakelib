# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
import os
import unittest
import logging
from nose.tools import raises
from snakemake.workflow import Workflow
from snakemake.exceptions import UnknownRuleException, NoRulesException
from snakemakelib.rules import create_rule_from_existing

logging.basicConfig(level=logging.DEBUG)

class TestRules(unittest.TestCase):
    """Test rules"""
    def setUp(self):
        self.workflow = Workflow("foo")
        name = self.workflow.add_rule(name="bar")
        self.workflow._rules["bar"].set_output(*(), **{"foo":"bar"})
        self.workflow._rules["bar"].set_params(*(), **{"cmd":"foo", "options":["foo", "bar"]})
    
    def test_create_rule(self):
        self.assertListEqual(["bar"], [x.name for x in self.workflow.rules])
        create_rule_from_existing(name="foo", template="bar", workflow=self.workflow)
        self.assertListEqual(["bar", "foo"], [x.name for x in self.workflow.rules])

    @raises(NoRulesException)
    def test_create_rule_empty_workflow(self):
        create_rule_from_existing(name="foo", template="bar", workflow=Workflow("foo"))

    @raises(AssertionError)
    def test_create_rule_wrong_workflow(self):
        create_rule_from_existing(name="foo", template="bar", workflow=None)

    @raises(UnknownRuleException)
    def test_create_rule_wrong_template(self):
        create_rule_from_existing(name="bar", template="foo", workflow=self.workflow)

    def test_create_rule_add_output(self):
        create_rule_from_existing(name="foo", template="bar", workflow=self.workflow, **{'output':((), {'bar':'foo'})})
        self.assertDictEqual({'bar':'foo'}, dict(self.workflow.get_rule("foo").output))

    def test_create_rule_add_params(self):
        create_rule_from_existing(name="foo", template="bar", workflow=self.workflow, **{'params':((), {'cmd':'bar', 'options':['foo']})})
        self.assertDictEqual({'cmd':'bar', 'options':['foo']}, dict(self.workflow.get_rule("foo").params))
