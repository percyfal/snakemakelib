# Copyright (C) 2014 by Per Unneberg
"""
Configuration module
"""
import os
from snakemakelib.log import LoggerManager
from collections import Mapping

smllogger = LoggerManager().getLogger(__name__)


def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(__file__)

def sml_rules_path():
    return os.path.join(os.path.dirname(__file__), "rules")

def sml_templates_path():
    return os.path.join(os.path.dirname(__file__), "_templates")



def update_config(d, u, overwrite_config=None):
    """Recursively update dictionary d with u.

    See
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for details.

    Args:
      d (dict): dictionary to update
      u (dict): dictionary whose items will overwrite those in d
      overwrite_config (dict): configuration passed via command line; keys in this dictionary will not be overwritten

    Returns:
      dict: updated dictionary

    """
    if overwrite_config is None:
        overwrite_config = dict()
    for (key, value) in u.items():
        if (isinstance(value, Mapping)):
            d[key]= update_config(d.get(key, {}), value, overwrite_config.get(key, {}))
        else:
            if not key in overwrite_config:
                d[key] = u[key]
            if isinstance(d[key], str):
                d[key] = os.path.expandvars(d[key])
    return d
