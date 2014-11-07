# Copyright (C) 2014 by Per Unneberg

global sml_config

sml_config = {'section1': 'value',
              'section2': {'subsection1':'value1',
                           'subsection2':'value2'}}

class BaseConfig(dict):
    _sections = []
    def __init__(self, *args, **kwargs):
        self._sections = [kk for k in args for kk in list(k)] + list(kwargs)
        self.update(*args, **kwargs)

    def __setitem__(self, key, val):
        if not key in self._sections:
            raise KeyError("key '" + key + "' not found in configuration dictionary")
        dict.__setitem__(self, key, val)
