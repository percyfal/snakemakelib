# Copyright (C) 2014 by Per Unneberg
"""
Configuration module
"""
import os
import sys
import yaml
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

class DeprecatedException(Exception):
    pass

def init_sml_config(cfg):
    raise DeprecatedException("init_sml_config is deprecated; snakemakelib now makes use of snakemakes internal configuration variable. See https://github.com/percyfal/snakemakelib/wiki for more information")

def get_sml_config(section=None):
    raise DeprecatedException("get_sml_config is deprecated; snakemakelib now makes use of snakemakes internal configuration variable. See https://github.com/percyfal/snakemakelib/wiki for more information")

# TODO: rename default    
def update_sml_config(config_default):
    raise DeprecatedException("update_sml_config is deprecated; snakemakelib now makes use of snakemakes internal configuration variable. See https://github.com/percyfal/snakemakelib/wiki for more information")

class BaseConfig(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        self._sections = []
        self.update(*args, **kwargs)

    def _inspect_sections(self):
        """Walk through configuration object to make sure subsections are BaseConfig classes, not dictionaries"""
        for k,v in self.items():
            if isinstance(v, dict):
                if not isinstance(v, BaseConfig):
                    smllogger.debug("Updating key {k} to <BaseConfig> class".format(k=k))
                    self[k] = BaseConfig(v)
                self[k]._inspect_sections()

    def __setitem__(self, key, val):
        # if not key in self._sections:
        #     raise KeyError("section '{key}' not found in configuration dictionary".format(key=key))
        if not key in self._sections:
            self.add_section(key)
        if isinstance(val, dict) and not isinstance(val, BaseConfig):
            val = BaseConfig(val)
        dict.__setitem__(self, key, val)

    def __getitem__(self, k):
        param = None
        if isinstance(k, tuple):
            key, param = k
        else:
            key = k
        val = dict.__getitem__(self, key)
        if str(type(val)) == "<class 'function'>":
            if not param is None:
                return val(param)
            return val()
        else:
            return val

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
            smllogger.error("Section {section} already present in configuration ".format(section=section))
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
            smllogger.error("Error: no such section {} in config; returning entire config".format(section))
            return self

    @property
    def sections(self):
        """Return sections"""
        return self._sections

def load_sml_config(config, cfg_file=None):
    """Load snakemakelib configuration file(s).

    Will search for configuration files in following order:
    
    1. ~/.smlconf.yaml - a personal site-wide configuration file
    2. ./smlconf.yaml - a standard configuration file residing in the
        same directory as the Snakefile
    3. cfg_file, if provided

    Args:
      config: snakemake configuration object
      cfg_file: custom configuration file to load

    Returns:
      config: updated configuration
    """
    config = BaseConfig(config)
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
        config.update(cfg)
    return config

def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(os.path.dirname(__file__))

def sml_rules_path():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "rules")

def sml_templates_path():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates")


##############################
# New configuration functions
# Work directly on global config object
##############################
def update_snakemake_config(config, update_config):
    """Update configuration object.

    Args:
        config: snakemake global configuration object
        update_config: configuration object of type <dict> or <BaseConfig>
    """    
    if not isinstance(config, BaseConfig):
        raise TypeError(
            """config object is not a <BaseConfig> object;
            you *must* do a 'config = load_sml_config(config)' statement prior to including any
            snakemakelib rules""")
    try:
        update_config = BaseConfig(update_config)
    except:
        raise TypeError("Failed to convert update_config object to <BaseConfig> object; recieved {cd}".format(cd=type(update_config)))
    config = _update_snakemake_config(config, update_config)
    return config

def _update_snakemake_config(config, update_config):
    """Update snakemake global configuration object. The default object is
    defined in the preamble of rules files and contains sensitive
    default settings.

    While traversing the config dictionary, make sure that the
    configuration settings correspond to those in the default
    dictionary, i.e. that if a section in default points to another
    BaseConfig object, then the config object should also point to a
    BaseConfig object.

    Loops through items in update_config and updates the config
    configuration. There are two cases:

    1. If the key/value pair is not present in config the value in
       update_config is used to set the coresponding value in
       config

    2. Else, use the set value in config.

    This procedure ensures that if a section is undefined in
    config, a default value will always be present.

    Args:
      update_config: configuration object to update global config with

    Returns:
      config: updated configuration

    """
    if not isinstance(update_config, BaseConfig):
        try:
            update_config = BaseConfig(update_config)
        except:
            raise TypeError("Failed to convert update_config to <BaseConfig> object; recieved {cd}".format(cd=type(update_config)))
    if config is None:
        return update_config
    if not type(config) == type(update_config):
        raise TypeError("config type {config}, default type {update_config}; configuration entry is not of type {type}".format(config=type(config), update_config=type(update_config), type=type(update_config)))
    # Loop sections
    for (section, value) in update_config.items():
        if not config.has_section(section):
            config.add_section(section)
        if (not isinstance(dict(update_config)[section], BaseConfig)):
            smllogger.debug("Key is not type BaseConfig: got type '{type}'".format(type = type(dict(update_config)[section])))
            # if config has no value set to default
            if config.get(section) is None:
                config[section] = dict(update_config)[section]
            # else make sure variable is expanded
            else:
                if isinstance(dict(config)[section], str):
                    smllogger.debug("expanding variables in config string '{s}', if present".format(s=dict(config)[section]))
                    config[section] = os.path.expandvars(config[section])
        else:
            config[section] = _update_snakemake_config(config[section], dict(update_config)[section])
    return config
