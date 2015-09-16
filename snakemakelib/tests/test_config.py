# Copyright (C) 2014 by Per Unneberg
# pylint: disable=R0904
import unittest
import logging
from nose.tools import raises
from snakemake.utils import update_config
import pytest

logging.basicConfig(level=logging.DEBUG)


def test_update_config():
    a = {'foo':{'bar':{'foo': 'bar', 'bar':'foo'}}}
    b = {'foo':{'bar':{'foo': 'barfoo'}}}
    update_config(a, b)
    c = a
    update_config(b, a)
    d = b
    assert c == {'foo': {'bar': {'foo': 'barfoo', 'bar': 'foo'}}}
    assert d == {'foo': {'bar': {'bar': 'foo', 'foo': 'barfoo'}}}

def test_update_config_with_string():
    with pytest.raises(AttributeError):
        update_config({}, "foo")

    


