# Copyright (c) 2014 Per Unneberg
import os
import re
import csv

# http://stackoverflow.com/questions/2170900/get-first-list-index-containing-sub-string-in-python
def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i
    return -1

def _raw(x):
    return (x, None)

def _convert_input(x):
    if re.match("^[0-9]+$", x):
        return int(x)
    elif re.match("^[0-9,.]+$", x):
        return float(x.replace(",", "."))
    else:
        return str(x)

def _read_picard_metrics(f):
    with open(f) as fh:
        data = fh.readlines()
        # Find histogram line
        i_hist = index_containing_substring(data, "## HISTOGRAM")
        if i_hist == -1:
            i = len(data)
        else:
            i = i_hist
        metrics = [[_convert_input(y) for y in x.rstrip("\n").split("\t")] for x in data[0:i] if not re.match("^[ #\n]", x)]
        if i_hist == -1:
            return (metrics, None)
        hist = [[_convert_input(y) for y in x.rstrip("\n").split("\t")] for x in data[i_hist:len(data)] if not re.match("^[ #\n]", x)]
    return (metrics, hist)

def _indent_texttable_for_rst(ttab, indent=4, add_spacing=True):
    """Texttable needs to be indented for rst.

    :param ttab: texttable object
    :param indent: indentation (should be 4 *spaces* for rst documents)
    :param add_spacing_row: add additional empty row below class directives

    :returns: reformatted texttable object as string
    """
    output = ttab.draw()
    new_output = []
    for row in output.split("\n"):
        new_output.append(" " * indent + row)
        if re.search('.. class::', row):
            new_row = [" " if x != "|" else x for x in row]
            new_output.append(" " * indent + "".join(new_row))
    return "\n".join(new_output)

# Convenience functions
def align2dict(align_metrics_tuple):
    """Get the alignment data in one big table.

    Args:
      list of (sample, AlignMetrics) tuples

    Returns:
      OrderedDict of align metrics
    """
    amcsv = []
    for t in align_metrics_tuple:
        amcsv = amcsv + t.metrics(as_csv=True)

    return amcsv
#df = [row for row in csv.DictReader(t.metrics(as_csv=True))]


# For now: extension maps to tuple (label, description). Label should
# be reused for analysis definitions
EXTENSIONS={'.align_metrics':('align', 'alignment', _read_picard_metrics),
            '.hs_metrics':('hs', 'hybrid selection', _read_picard_metrics),
            '.dup_metrics':('dup', 'duplication metrics', _read_picard_metrics),
            '.insert_metrics':('insert', 'insert size', _read_picard_metrics),
            '.eval_metrics':('eval', 'snp evaluation', _raw)
            }

class PicardMetrics(object):
    """PicardMetrics: class for reading/storing metrics from one picard
    metrics file.

    Args:
      pmid (str): unique identifier
      file (str): metrics file to read

    """
    def __init__(self, pmid, f):
        self._metrics = None
        self._hist = None
        self._file = f
        self._pmid = pmid
        (_, self._metrics_type) = (os.path.splitext(f))
        self._read_metrics()

    def _read_metrics(self):
        """Read metrics"""
        (self._metrics, self._hist) = EXTENSIONS[self._metrics_type][2](self._file)

    def metrics(self, as_csv=False):
        if as_csv:
            return [",".join([str(y) for y in x]) for x in self._metrics]
        return self._metrics

    def hist(self, as_csv=False):
        if as_csv:
            return [",".join([str(y) for y in x]) for x in self._hist]
        return self._hist

    @property
    def id(self):
        return self._pmid

    @property
    def type(self):
        return self._metrics_type

    @property
    def file(self):
        return self._file

class AlignMetrics(PicardMetrics):
    def __init__(self, pmid, f):
        PicardMetrics.__init__(self, pmid, f)

class InsertMetrics(PicardMetrics):
    def __init__(self, pmid, f):
        PicardMetrics.__init__(self, pmid, f)

class DuplicationMetrics(PicardMetrics):
    def __init__(self, pmid, f):
        PicardMetrics.__init__(self, pmid, f)
        
class HsMetrics(PicardMetrics):
    def __init__(self, pmid, f):
        PicardMetrics.__init__(self, pmid, f)
    

class PicardMetricsCollection(object):
    """Collection of PicardMetrics objects. NB: name currently misleading
    as it is *not* an iterable!
    
    Args:

      mlist (list): list of tuples consisting of two elements (sample, metrics_file)

    """
    def __init__(self, mlist):
        self._metrics = {}
        self._mlist = mlist
        self._collect_metrics()

    def _collect_metrics(self):
        for (sid, fn) in self._mlist:
            pm = PicardMetrics(sid, fn)
            # FIXME: could be more of same type, e.g. several
            # dup_metrics before and after duplicate removal
            self._metrics[pm.type] = pm
            
    def metrics(self, as_csv=False):
        return {pm.type : pm.metrics(as_csv=as_csv) for pm in sorted(self._metrics.values(), key=lambda x: x.type)}

    def hist(self, as_csv=False):
        return {pm.type : pm.hist(as_csv=as_csv) for pm in sorted(self._metrics.values(), key=lambda x: x.type)}

    def idlist(self):
        return [pm.id for pm in sorted(self._metrics.values(), key=lambda x: x.type)]

        
