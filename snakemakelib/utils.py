# Copyright (c) 2014 Per Unneberg
import os
import shutil
from datetime import datetime, date
from snakemakelib.log import LoggerManager
from snakemakelib.stat import is_installed
from jinja2 import Environment, PackageLoader

smllogger = LoggerManager().getLogger(__name__)

# TODO: move this elsewhere
SmlTemplateEnv = Environment(loader = PackageLoader("snakemakelib", "_templates"))
SmlTemplateEnv.globals.update(zip=zip)

def utc_time():
    """Make an utc_time with appended 'Z'"""
    return str(datetime.utcnow()) + 'Z'

def isoformat(s=None):
    """Return isoformat date from string"""
    if s is None:
        return
    # Assume YYMMDD format
    if len(s) == 6:
        (YY, MM, DD) = (s[0:2], s[2:4], s[4:6])
        return date(int("20{YY}".format(YY=YY)), int(MM.lstrip("0")), int(DD)).isoformat()

# http://stackoverflow.com/questions/2556108/how-to-replace-the-last-occurence-of-an-expression-in-a-string
def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

## From bcbb
def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not os.path.exists(dname):
        try:
            os.makedirs(dname)
        except OSError:
            if not os.path.isdir(dname):
                raise
    else:
        smllogger.warning("Directory {} already exists; not making directory".format(dname))
    return dname

class NotInstalledError(Exception):
    """Error thrown if program/command/application cannot be found in path

    Args:
      msg (str): String described by exception
      code (int, optional): Error code, defaults to 2.
    
    """
    def __init__(self, msg, code=2):
        self.msg = msg
        self.code = code

def set_cmd(home, cmd, module):
    """Set the command, checking if the program is installed in the
    process.

    Args:
      home (str): path to application
      cmd (str): the actual command name
      module (str): the calling module

    Returns:
      str: full path to the command

    Example:

        cmd = set_cmd("/path/to/cmd", "helloworld")
        print(cmd)
        # prints /path/to/cmd/helloworld

    Raises:
      NotInstalledError: if program not found in path

    """
    if home:
        cmd = os.path.join(home, cmd)
    else:
        try:
            cmd = shutil.which(cmd)
        except:
            pass

    if not is_installed(cmd):
        raise NotInstalledError("\n{module}: {prog} not installed or not in PATH\n".format(module=module, prog=cmd))

    return cmd
