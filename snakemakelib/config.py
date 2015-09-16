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

class DeprecatedException(Exception):
    pass

def expandvars_in_config(d):
    """Expand variables in config

    Traverse config and use os.path.expandvars to expand strings.

    """
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

    DEPRECATED!
    """
    raise DeprecatedException("""snakemakelib.config.load_sml_config has been deprecated: use the 'configfile: "configfile"' syntax instead""")

def update_config(config, overwrite_config):
    """Recursively update dictionary config with overwrite_config.
    
    Moved to snakemake.utils.
    """
    raise DeprecatedException("""snakemakelib.config.update_config has been deprecated as it is now part of snakemake: use snakemake.utils.update_config""")
