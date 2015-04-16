# Copyright (c) 2014 Per Unneberg
from datetime import datetime, date

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
