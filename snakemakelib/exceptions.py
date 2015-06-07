# Copyright (C) 2015 by Per Unneberg


class NotInstalledError(Exception):
    """Error thrown if program/command/application cannot be found in path

    Args:
      msg (str): String described by exception
      code (int, optional): Error code, defaults to 2.

    """
    def __init__(self, msg, code=2):
        self.msg = msg
        self.code = code


class SamplesException(Exception):
    """Error thrown if samples missing or wrong number.

    Args:
      msg (str): String described by exception
      code (int, optional): Error code, defaults to 2.

    """
    def __init__(self, msg, code=2):
        self.msg = msg
        self.code = code


class OutputFilesException(Exception):
    """Error thrown if outputfiles missing or wrong number.

    Args:
      msg (str): String described by exception
      code (int, optional): Error code, defaults to 2.

    """
    def __init__(self, msg, code=2):
        self.msg = msg
        self.code = code
