# Copyright (c) 2014 Per Unneberg
import logging
import logging.config
import yaml
import os

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)

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
    _ch = logging.StreamHandler()
    _formatter = logging.Formatter(_fmt)
    _ch.setFormatter(_formatter)
    _has_loaded_config = False

    def __init__(self, *args, **kwargs):
        if not LoggerManager._has_loaded_config:
            # Add snakemakelib root handler
            LoggerManager._loggers['snakemakelib'] = logging.getLogger('snakemakelib')
            LoggerManager._loggers['snakemakelib'].setLevel(logging.WARNING)
            LoggerManager._loggers['snakemakelib'].addHandler(LoggerManager._ch)
            self._load_config()
            LoggerManager._has_loaded_config = True

    def _load_config(self):
        conf = {}
        if os.path.exists("logconf.yaml"):
            with open ("logconf.yaml", "r") as fh:
                conf = yaml.load(fh)
        else:
            conf = self._load_sml_config()
        if conf:
            logging.config.dictConfig(conf)
     
    def _load_sml_config(self):
        cfg = {}
        for fn in [os.path.join(os.getenv("HOME"), ".smlconf.yaml"),
                   os.path.join(os.curdir, "smlconf.yaml")]:
            if (fn is None):
                continue
            if not os.path.exists(fn):
                continue
            with open(fn, "r") as fh:
                cfg_tmp = yaml.load(fh)
            if cfg_tmp:
                cfg.update(cfg_tmp)
        return cfg.get("logging", {})

    @staticmethod
    def getLogger(name=None):
        if not name:
            smllogger = logging.getLogger()
            return logging.getLogger()
        elif name not in LoggerManager._loggers.keys():
            LoggerManager._loggers[name] = logging.getLogger(str(name))
            LoggerManager._loggers[name].addHandler(LoggerManager._ch)
            return LoggerManager._loggers[name]
        else:
            logging.warn("Trying to get already existing logger")
