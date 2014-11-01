# Copyright (C) 2014 by Per Unneberg

"""
Utility functions for snakemakelib.
"""
import os

def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(os.path.dirname(__file__))

def sml_rules_path():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "rules")
