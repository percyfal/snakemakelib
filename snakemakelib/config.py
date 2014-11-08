# Copyright (C) 2014 by Per Unneberg
"""
Configuration module
"""

import os
from snakemake.logging import logger

def update_config(config, d):
    for (section, value) in d.items():
        config[section] = config.get(section, d[section])
        if (not isinstance(d[section], dict)):
            # config key has value - can be modified via commandline
            config[section] = config.get(section, value)
        else:
            for (subsection, option) in d[section].items():
                if (isinstance(option, dict)):
                    config[section][subsection] = config[section].get(subsection, {})
                    for (param, value) in option.items():
                        config[section][subsection][param] = config[section][subsection].get(param, value)
                else:
                    config[section][subsection] = config[section].get(subsection, option)
    return config

def update_sml_config(config, d):
    """Update snakemakelib configuration object

    :param config: configuration to update
    :param d: config with new values

    :return: updated configuration object
    """
    for (section, value) in d.items():
        config[section] = config.get(section, d[section])
        if (not isinstance(d[section], dict)):
            # config key has value - can be modified via commandline
            config[section] = config.get(section, value)
        else:
            for (subsection, option) in d[section].items():
                if (isinstance(option, dict)):
                    config[section][subsection] = config[section].get(subsection, {})
                    for (param, value) in option.items():
                        config[section][subsection][param] = config[section][subsection].get(param, value)
                else:
                    config[section][subsection] = config[section].get(subsection, option)
    return config
    


def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(os.path.dirname(__file__))

def sml_rules_path():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "rules")
                            

global sml_config

sml_config = {'section1': 'value',
              'section2': {'subsection1':'value1',
                           'subsection2':'value2'}}

class BaseConfig(dict):
    _sections = []

    def _inspect_sections(self):
        """Walk through configuration object to make sure subsections are BaseConfig classes, not dictionaries"""
        for k,v in self.items():
            logger.info("{k} : {v}".format(k=k, v=v))
            if isinstance(v, dict):
                if not isinstance(v, BaseConfig):
                    raise TypeError("section {k}; dictionary {v} must be instance of <BaseConfig> class".format(k=k,v=v) )
                v._inspect_sections()

    def __init__(self, *args, **kwargs):
        self._sections = [kk for k in args for kk in list(k)] + list(kwargs)
        self.update(*args, **kwargs)
        # Make sure all dict subsections are BaseConfig objects
        self._inspect_sections()

    def __setitem__(self, key, val):
        if not key in self._sections:
            raise KeyError("section '{key}' not found in configuration dictionary".format(key=key))
        if isinstance(val, dict) and not isinstance(val, BaseConfig):
            raise TypeError("dictionary {val} must be instance of <BaseConfig> class".format(val=val) )
        dict.__setitem__(self, key, val)

    def add_section(self, section):
        """Add new section to configuration object."""
        if not (isinstance(section, str)):
            raise TypeError("argument 'section' must be of type <str>")
        if section in self._sections:
            logger.error("Section {section} already present in configuration ".format(section=section))
            return
        self._sections.append(section)
        self[section] = None

    @property
    def sections(self):
        """Return sections"""
        return self._sections
