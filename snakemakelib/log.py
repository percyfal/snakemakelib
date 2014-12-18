# Copyright (c) 2014 Per Unneberg
import logging

def setup_logging():
    # use a variable in the function object to determine if it has run before
    if getattr(setup_logging, "has_run", False):
        return

    logger = logging.getLogger('snakemakelib')
    logger.setLevel(logging.DEBUG)

    streamHandler = logging.StreamHandler()
    streamHandler.setLevel(logging.DEBUG)

    file_format = "%(asctime)s (%(levelname)s) %(name)s : " + \
        "%(message)s"

    formatter = logging.Formatter(file_format)
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)
    setup_logging.has_run = True

def get_logger():
    logger = logging.getLogger('snakemakelib')
    if not logger.handlers:
        setup_logging()
    return logger
