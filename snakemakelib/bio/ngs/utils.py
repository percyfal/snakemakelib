# Copyright (C) 2015 by Per Unneberg
import re
import os
from snakemakelib.bio.ngs.regexp import RegexpDict

def find_files(regexp, path = os.curdir, search=False):
    """Find files in path that comply with a regular expression.

    Args:
      regexp: regular expression object of class <RegexpDict> or string
      path:   path to search
      search: use re.search instead of re.match for pattern matching

    Returns:
      flist: list of file names, prepended with root path
    """
    if isinstance(regexp, RegexpDict):
        r = regexp.re
    else:
        r = re.compile(regexp)
    re_fn = r.search if search else r.match
    flist = []
    for root, dirs, files in os.walk(path):
        flist += [os.path.join(root, x) for x in files if re_fn(x)]
    return sorted(flist)
