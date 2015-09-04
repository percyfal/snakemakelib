# Copyright (C) 2014 by Per Unneberg
"""
Configuration module
"""
import os
import yaml
from snakemakelib.log import LoggerManager
from collections import Mapping

smllogger = LoggerManager().getLogger(__name__)

SNAKEMAKELIB_PATH = os.path.dirname(__file__)
SNAKEMAKELIB_RULES_PATH = os.path.join(SNAKEMAKELIB_PATH, "rules")
SNAKEMAKELIB_TEMPLATES_PATH = os.path.join(SNAKEMAKELIB_PATH, "_templates")

def _expandvars(d):
    def _update(d):
        for k, v in d.items():
            if isinstance(v, Mapping):
                r = _update(d.get(k, {}))
                d[k] = r
            else:
                if isinstance(d[k], str):
                    d[k] = os.path.expandvars(d[k])
        return d
    _update(d)

def load_sml_config(config, config_file=None, expandvars=True):
    """Load snakemakelib configuration file(s).
    Will search for configuration files in following order:
    
    1. ~/.smlconf.yaml - a personal site-wide configuration file
    2. ./smlconf.yaml - a standard configuration file residing in the
        same directory as the Snakefile
    3. config_file, if provided

    Args:

      config (dict): snakemake configuration object
      config_file (str): custom configuration file to load
      expandvars (bool): do automatic expansion of environment variables

    """
    for fn in [os.path.join(os.getenv("HOME"), ".smlconf.yaml"),
               os.path.join(os.curdir, "smlconf.yaml"),
               config_file]:
        if (fn is None):
            continue
        if not os.path.exists(fn):
            continue
        smllogger.info("Loading configuration from {}".format(fn))
        with open(fn, "r") as fh:
            overwrite_config = yaml.load(fh)
        smllogger.info("Read configuration from {}".format(fn))
        update_config(config, overwrite_config)
    if expandvars:
        _expandvars(config)

def update_config(config, overwrite_config):
    """Recursively update dictionary config with overwrite_config.

    See
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for details.

    Args:
      config (dict): dictionary to update
      overwrite_config (dict): dictionary whose items will overwrite those in config

    """
    def _update(d, u):
        for (key, value) in u.items():
            if (isinstance(value, Mapping)):
                d[key]= _update(d.get(key, {}), value)
            else:
                d[key] = u[key]
        return d
    config = _update(config, overwrite_config)
