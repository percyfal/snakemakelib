# Copyright (C) 2015 by Per Unneberg
import os
import shutil
import re

def is_compressed(f, comp_re):
    """Ascertain whether a file is compressed or not.

    Args:
      f: filename
      comp_re: regular expression of compression suffixes

    returns:
      boolean
    """
    return not re.search(comp_re, f) is None

def is_installed(prog):
    if not shutil.which(prog) is None:
        return True
    if os.path.exists(prog):
        return True
    return False
