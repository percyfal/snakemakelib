# Copyright (C) 2014 by Per Unneberg
"""
Configuration module
"""
import os
from snakemake.logging import logger

class BaseConfig(dict):
    def _inspect_sections(self):
        """Walk through configuration object to make sure subsections are BaseConfig classes, not dictionaries"""
        for k,v in self.items():
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

    def get_section(self, section):
        """Return section from config"""
        if not (isinstance(section, str)):
            raise TypeError("argument 'section' must be of type <str>")
        try:
            return self[section]
        except KeyError:
            logger.error("Error: no such section {} in config; returning entire config".format(section))
            return self

    @property
    def sections(self):
        """Return sections"""
        return self._sections

# Global configuration variable
__sml_config__ = BaseConfig({})

def init_sml_config(cfg):
    """Initialize sml configuration.

    Args:
        cfg: A configuration object of type <BaseConfig>
    """
    global __sml_config__
    __sml_config__ = BaseConfig({})
    __sml_config__ =  _update_sml_config(__sml_config__, cfg)

def get_sml_config(section=None):
    if not section is None:
        return __sml_config__.get_section(section)
    return __sml_config__

# TODO: rename default    
def update_sml_config(default):
    """Update sml configuration object.

    Loops through items in default configuration and updates the cfg
    configuration. There are two cases:

    1. If the key/value pair is not present in custom_cfg the value in
       default is used to set the coresponding value in custom_cfg
    2. Else, use the set value in custom_cfg.

    This procedure ensures that if a section is undefined in
    custom_cfg, a defualt value will always be present.

    Args:
        custom_cfg: A configuration object with custom settings of
                    type <BaseConfig>
        default: default configuration object of type <BaseConfig>
    """    
    global __sml_config__
    _update_sml_config(__sml_config__, default)

def _update_sml_config(config, default):
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
    if not isinstance(default, BaseConfig):
        raise TypeError("default config must be instance <class 'snakemakelib.config.BaseConfig'>; found {type}".format(type=type(default)) )        
    if not type(config) == type(default):
        raise TypeError("config type {config}, default type {default}; configuration entry is not of type {type}".format(config=type(config), default=type(default), type=type(default)))
    # Loop sections
    for (section, value) in default.items():
        if not config.has_section(section):
            config.add_section(section)
        if (not isinstance(default[section], BaseConfig)):
            # if config has no value set to default
            if config.get(section) is None:
                config[section] = default[section]
        else:
            config[section] = _update_sml_config(config[section], default[section])
    # if not set(list(config)).issuperset(set(list(default))):
    #     print(list(config))
    #     print(list(default))
    #     # TODO: make this a warning
    #     logger.info("Sections {} not defined in default rules configuration; this setting will not affect rules behaviour".format(list(set(list(config)).difference(set(list(default))))))
    return config

def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(os.path.dirname(__file__))

def sml_rules_path():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "rules")

