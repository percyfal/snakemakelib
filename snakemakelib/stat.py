# Copyright (C) 2015 by Per Unneberg
import os
import shutil
import re
from snakemakelib.config import get_sml_config

def is_compressed(f):
    """Ascertain whether a file is compressed or not"""
    cfg = get_sml_config('comp.settings')
    return not re.search(cfg['compression']['re'], f) is None

def is_installed(prog):
    if not shutil.which(prog) is None:
        return True
    if os.path.exists(prog):
        return True
    return False
