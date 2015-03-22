# Copyright (c) 2014 Per Unneberg
import logging
import logging.config
import yaml
import os
import snakemakelib.config

# See http://stackoverflow.com/questions/15727420/using-python-logging-in-multiple-modules
class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances.keys():
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class LoggerManager(object):
    __metaclass__ = Singleton

    _loggers = {}
    _fmt = "%(asctime)s (%(levelname)s) %(name)s :  %(message)s"

    def __init__(self, *args, **kwargs):
        if os.path.exists("logconf.yaml"):
            with open ("logconf.yaml", "r") as fh:
                conf = yaml.load(fh)
            logging.config.dictConfig(conf)
        else:
            snakemakelib.config.load_sml_config()
            conf = snakemakelib.config.get_sml_config('logging')
            logging.config.dictConfig(conf)
        LoggerManager._loggers['snakemakelib'] = logging.getLogger('snakemakelib').addHandler(logging.NullHandler())


    @staticmethod
    def getLogger(name=None):
        if not name:
            logging.basicConfig(format=LoggerManager._fmt, level=logging.INFO)
            return logging.getLogger()
        elif name not in LoggerManager._loggers.keys():
            logging.basicConfig(format=LoggerManager._fmt, level=logging.INFO)
            LoggerManager._loggers[name] = logging.getLogger(str(name))
        return LoggerManager._loggers[name]    
