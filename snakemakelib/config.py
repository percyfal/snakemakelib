# Copyright (C) 2014 by Per Unneberg
"""
Configuration module
"""
import os
import yaml
from snakemakelib.log import LoggerManager
from collections import Mapping

smllogger = LoggerManager().getLogger(__name__)

def load_sml_config(config, cfg_file=None):
    """Load snakemakelib configuration file(s).
    Will search for configuration files in following order:
    
    1. ~/.smlconf.yaml - a personal site-wide configuration file
    2. ./smlconf.yaml - a standard configuration file residing in the
        same directory as the Snakefile
    3. cfg_file, if provided
    Args:
      config (dict): snakemake configuration object
      cfg_file (str): custom configuration file to load
    Returns:
      config (dict): updated configuration
    """
    for fn in [os.path.join(os.getenv("HOME"), ".smlconf.yaml"),
               os.path.join(os.curdir, "smlconf.yaml"),
               cfg_file]:
        if (fn is None):
            continue
        if not os.path.exists(fn):
            continue
        smllogger.info("Loading configuration from {}".format(fn))
        with open(fn, "r") as fh:
            cfg = yaml.load(fh)
        smllogger.info("Read configuration from {}".format(fn))
        config = update_config(config, cfg)
    return config


def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(__file__)

def sml_rules_path():
    return os.path.join(os.path.dirname(__file__), "rules")

def sml_templates_path():
    return os.path.join(os.path.dirname(__file__), "_templates")

def update_config(d, u):
    """Recursively update dictionary d with u.

    See
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for details.

    Args:
      d (dict): dictionary to update
      u (dict): dictionary whose items will overwrite those in d

    Returns:
      dict: updated dictionary

    """
    for (key, value) in u.items():
        if (isinstance(value, Mapping)):
            d[key]= update_config(d.get(key, {}), value)
        else:
            d[key] = u[key]
            if isinstance(d[key], str):
                d[key] = os.path.expandvars(d[key])
    return d
