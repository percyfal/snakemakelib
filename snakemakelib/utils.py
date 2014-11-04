# Copyright (C) 2014 by Per Unneberg

"""
Utility functions for snakemakelib.
"""
import os
import collections

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



def sml_path():
    return os.path.dirname(__file__)

def sml_base_path():
    return os.path.dirname(os.path.dirname(__file__))

def sml_rules_path():
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "rules")
                            
