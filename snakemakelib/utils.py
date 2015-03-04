# Copyright (c) 2014 Per Unneberg
import re
from datetime import datetime
from snakemakelib.config import get_sml_config

def utc_time():
    """Make an utc_time with appended 'Z'"""
    return str(datetime.utcnow()) + 'Z'

def is_compressed(f):
    """Ascertain whether a file is compressed or not"""
    cfg = get_sml_config('comp.settings')
    return not re.search(cfg['compression']['re'], f) is None
