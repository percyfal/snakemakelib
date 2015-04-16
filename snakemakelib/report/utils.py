# Copyright (c) 2014 Per Unneberg
import re
import string
import math
import itertools
from datetime import datetime
from collections import OrderedDict

class Template(string.Formatter):
    _suffix = {'-3':('n', 10**(-9)), '-2':('u', 10**(-6)), '-1':('m', 10**(-3)), '0':('', 1), '1':('k', 10**3), '2':('M', 10**6), '3':('G', 10**9), '4':('T', 10**12), '5':('P', 10**15)}
    def format_field(self, value, spec):
        sfx = ""
        if spec.endswith('h'):
            if not value == 0:
                spec = spec[:-1] + 'f'
                n = (math.floor(math.log(value,10) / 3))
                value = value / self._suffix[str(n)][1]
                sfx = self._suffix[str(n)][0]
            else:
                spec = 'd'
        return super(Template, self).format_field(value, spec) + sfx

def recast(x, strpfmt="%b %d %H:%M:%S"):
    x = x.rstrip().lstrip()
    if re.match('^[0-9]+$', x):
        return int(x)
    elif re.match('^[0-9]+[,\.][0-9]+$', x):
        return float(x.replace(",", "."))
    elif re.search("%", x):
        return float(x.replace(",", ".").replace("%", ""))
    else:
        try:
            dateobj = datetime.strptime(x, strpfmt)
            return dateobj
        except:
            return x

# Replace whitespace with underscore, convert percent characters to PCT
def trim_header(x, underscore=False, percent=False):
    return x.lstrip().rstrip().replace(" ", "_" if underscore else " ").replace("%", "PCT" if percent else "%").replace(",", "_" if underscore else " ")


def group_samples(samples, grouping="sample"):
    """Group samples by sample or sample run.

    Args:
      samples: list of namedtuples with fields (sample_id, project_id)
      grouping: what to group by

    Returns:
      groups: dictionary of grouped items
    """
    groups = OrderedDict()
    if grouping == "sample":
        for k,g in itertools.groupby(samples, key=lambda x:x.sample_id()):
            groups[k] = list(g)
    return groups
