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



def update_sml_config(config, default):
    """Update snakemakelib configuration object. The default object is
defined in the preamble of rules files and contains sensitive default
settings. 

    While traversing the config dictionary, make sure that the configuration settings correspond to those in the default dictionary, i.e. that if a section in default points to another BaseConfig object, then the config object should also point to a BaseConfig object.

    :param config: configuration to update
    :param default: config with default values

    :return: updated configuration object
    """
    if config is None:
        return default
    if not type(config) == type(default):
        raise TypeError("config {config}, default {default}; configuration entry is not of type {type}".format(config=config, default=default, type=type(default)))
    if isinstance(default, BaseConfig):
        if not set(list(config)).issuperset(set(list(default))):
            # TODO: make this a warning
            logger.info("Unique keys in config: {}".format(list(set(list(config)).difference(set(list(default))))))

    # Loop sections
    for (section, value) in default.items():
        if not config.has_section(section):
            config.add_section(section)
        if (not isinstance(default[section], BaseConfig)):
            # if config has no value set to default
            if config.get(section) is None:
                config[section] = default[section]
        else:
            config[section] = update_sml_config(config[section], default[section])
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
    def _inspect_sections(self):
        """Walk through configuration object to make sure subsections are BaseConfig classes, not dictionaries"""
        for k,v in self.items():
            logger.info("{k} : {v}".format(k=k, v=v))
            if isinstance(v, dict):
                if not isinstance(v, BaseConfig):
                    raise TypeError("section {k}; dictionary {v} must be instance of <BaseConfig> class".format(k=k,v=v) )
                v._inspect_sections()

    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        self._sections = []
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        if not key in self._sections:
            raise KeyError("section '{key}' not found in configuration dictionary".format(key=key))
        if isinstance(val, dict) and not isinstance(val, BaseConfig):
            raise TypeError("dictionary {val} must be instance of <BaseConfig> class".format(val=val) )
        dict.__setitem__(self, key, val)

    def update(self, *args, **kwargs):
        self._sections += [kk for k in args for kk in list(k)] + list(kwargs)
        dict.update(self, *args, **kwargs)
        # Make sure all dict subsections are BaseConfig objects
        self._inspect_sections()

    def add_section(self, section):
        """Add new section to configuration object."""
        if not (isinstance(section, str)):
            raise TypeError("argument 'section' must be of type <str>")
        if section in self._sections:
            logger.error("Section {section} already present in configuration ".format(section=section))
            return
        self._sections.append(section)
        self[section] = None

    def has_section(self, section):
        """Check if section is already defined"""
        return section in self.sections

    @property
    def sections(self):
        """Return sections"""
        return self._sections
