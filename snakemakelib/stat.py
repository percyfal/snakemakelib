# Copyright (C) 2015 by Per Unneberg
import re
from snakemakelib.config import get_sml_config

def is_compressed(f):
    """Ascertain whether a file is compressed or not"""
    cfg = get_sml_config('comp.settings')
    return not re.search(cfg['compression']['re'], f) is None

