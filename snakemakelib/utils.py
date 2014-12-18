# Copyright (c) 2014 Per Unneberg
from datetime import datetime

def utc_time():
    """Make an utc_time with appended 'Z'"""
    return str(datetime.utcnow()) + 'Z'
