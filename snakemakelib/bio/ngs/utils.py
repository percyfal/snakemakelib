# Copyright (C) 2015 by Per Unneberg
import re
import os
from snakemakelib.bio.ngs.regexp import RegexpDict


def find_files(regexp, path=os.curdir, search=False, limit=None):
    """Find files in path that comply with a regular expression.

    Args:
      regexp (RegexpDict | str): regular expression object of class
                               <RegexpDict> or <str>
      path (str):   path to search
      search (bool): use re.search instead of re.match for pattern matching
      limit (dict): dictionary where keys correspond to regular expression
             grouping labels and values are lists that limit the
             returned pattern

    Returns:
      flist: list of file names, prepended with root path
    """
    if isinstance(regexp, RegexpDict):
        r = regexp.re
    else:
        if not regexp:
            return []
        r = re.compile(regexp)
    re_fn = r.search if search else r.match
    flist = []
    for root, dirs, files in os.walk(path):
        for x in files:
            m = re_fn(x)
            if m is None:
                continue
            if limit:
                if any([m.group(k) in limit[k] for k in limit.keys()]):
                    flist += [os.path.join(root, x)]
            else:
                flist += [os.path.join(root, x)]
    return sorted(flist)

def dict_to_R(d, as_string=True):
    """Translate a python dictionary to an R option string

    Args:
      d: dictionary

    Returns:
      A string representing an option string of the dictionary in R
    """
    outlist = []
    for (k, v) in d.items():
        if isinstance(v, list):
            if isinstance(v[0], int):
                outlist.append("{k} = c(".format(k=k) + ",".join("{vv}".format(vv=vv) for vv in v) + ")")
            else:
                outlist.append("{k} = c(".format(k=k) + ",".join("'{vv}'".format(vv=vv) for vv in v) + ")")     
        elif isinstance(v, dict):
            outlist.append("{k} = c(".format(k=k) + ",".join("{kk}={vv}".format(kk=kk, vv=vv) for (kk,vv) in v.items()) + ")")
        elif v is True:
            outlist.append("{k}=TRUE".format(k=k))
        elif v is False:
            outlist.append("{k}=FALSE".format(k=k))
        elif v is None:
            outlist.append("{k}=NULL".format(k=k))
        elif isinstance(v, int):
            outlist.append("{k}={v}".format(k=k, v=v))
        elif isinstance(v, float):
            outlist.append("{k}={v}".format(k=k, v=v))
        else:
            outlist.append("{k}='{v}'".format(k=k, v=v))
    return ",".join(outlist)

def dict_to_Rdict(d, as_string=True):
    """Translate a python dictionary to a dict with R-compatible entries

    Args:
      d: dictionary

    Returns:
      A dictionare where values comply with R
    """
    dout = {}
    for (k, v) in d.items():
        if isinstance(v, list):
            if isinstance(v[0], int):
                dout[k] = "c(" + ",".join("{vv}".format(vv=vv) for vv in v) + ")"
            else:
                dout[k] = "c(" + ",".join("'{vv}'".format(vv=vv) for vv in v) + ")"
        elif isinstance(v, dict):
            dout[k] = "c(" + ",".join("{kk}={vv}".format(kk=kk, vv=vv) for (kk,vv) in v.items()) + ")"
        elif v is True:
            dout[k] = "TRUE"
        elif v is False:
            dout[k] = "FALSE"
        elif v is None:
            dout[k] = "NULL"
        elif isinstance(v, int):
            dout[k] = v
        elif isinstance(v, float):
            dout[k] = v
        else:
            dout[k] = '{v}'.format(v=v)
    return dout
