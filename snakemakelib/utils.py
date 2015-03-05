# Copyright (c) 2014 Per Unneberg
import re
from datetime import datetime, date
from snakemakelib.config import get_sml_config

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

def is_compressed(f):
    """Ascertain whether a file is compressed or not"""
    cfg = get_sml_config('comp.settings')
    return not re.search(cfg['compression']['re'], f) is None
