# Copyright (c) 2014 Per Unneberg
import os
import sys
import re
import csv
import collections
from snakemakelib.report.utils import Template

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

def _make_unique(l):
    cnt = collections.Counter(l)
    luniq = []
    d = {}
    for c in l:
        if not c in d.keys():
            d[c] = 0
        else:
            d[c] += 1
        luniq.append("{c}.{sfx}".format(c=c, sfx=d[c]) if d[c]>0 else c)
    return luniq

class DataFrame(object):
    """Light weight data frame object

    A data frame is represented as an OrderedDict.


    Args:
      *args: if provided, must be a list of lists that upon initialization is converted to an OrderedDict. The first list is treated as a header.
    """
    def __init__(self, *args):
        if (len(args[0]) == 1):
            self._colnames = args[0]
            self._data = [collections.OrderedDict([(args[0][0], row[0])]) for row in args[1:]]
        else:
            reader = csv.DictReader([",".join([str(y) for y in x]) for x in args])
            self._colnames = reader.fieldnames
            self._data = [collections.OrderedDict([(k, row[k]) for k in self._colnames]) for row in reader]

    def __str__(self):
        return "{cls} object with {rows} rows, {columns} columns".format(cls=self.__class__, rows=self.dim[0], columns=self.dim[1])

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if len(self._data) > self.index:
            self.index += 1
            return self._data[self.index - 1]
        else:
            raise StopIteration

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._data]
        return self.__class__(*a)

    def __setitem__(self, key, val):
        _ = [row.update([(key, val)]) for row in self._data]

    def x(self, column=None, indices=None):
        column = self.colnames[0] if not column else column
        if indices:
            return [self._data[i][column] for i in indices]
        return [row[column] for row in self._data]
    
    # Copy definition of x
    y=x

    @property
    def colnames(self):
        return self._colnames

    @property
    def data(self):
        return self._data

    @property
    def dim(self):
        return (len(self._data), len(self._data[0]))
    
class PicardMetrics(object):
    """Generic class to store metrics section from Picard Metrics reports.
    See also class PicardHistMetrics for reports that provide metrics
    and histogram information.

    Args:
      *args: if provided, must be a list of lists that upon initialization is converted to an OrderedDict. The first list is treated as a header.
      filename (str): file from which to collect metrics
      identifier (str): unique identifier for object, e.g. sample or unit name

    Returns:
      An instance of class PicardMetrics
    """
    _format = collections.OrderedDict()
    _tp = Template()

    def __init__(self, *args, filename=None, identifier=None):
        self._id = identifier
        if filename is None and not args:
            raise ValueError("please supply either filename or args to instantiate class")
        self._set_vars(identifier, filename)
        if not self.filename is None and not args:
            (args, _) = _read_picard_metrics(self.filename)
        self._metrics = DataFrame(*args)

    def _set_vars(self, identifier, filename):
        self._id = identifier if not identifier is None else filename
        self._filename = str(filename)

    def __str__(self):
        return "{cls} object with a metrics field with {rows} rows, {columns} columns".format(cls=self.__class__, rows=self.metrics.dim[0], columns=self.metrics.dim[1])

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics._data]
        return self.__class__(*a, identifier=self.id, filename=self.filename)

    @property
    def metrics(self):
        return self._metrics
    
    @property 
    def id(self):
        return self._id

    @property
    def filename(self):
        return self._filename

    def _format_field(self, value, spec, ctype):
        if value == '?':
            spec = 's'
            ctype = str
        elif value == "":
            spec = 's'
            ctype = str
        value = ctype(value)
        return self._tp.format_field(value, spec)

    def summary(self, fmt = None, ctype = None, sep="\t", raw=False):
        columns = self.metrics.colnames
        fmt = {k:v[0] for (k,v) in list(self._format.items()) if k in columns} if fmt is None else {k:v for (k,v) in zip(columns, fmt)}
        ctype = {k:v[1] for (k,v) in list(self._format.items()) if k in columns} if ctype is None else {k:v for (k,v) in zip(columns, ctype)}
        if not fmt:
            raise ValueError("No format defined for {cls}; please use derived subclass".format(cls=__class__))
        if raw:
            return "\n".join([sep.join([x for x in columns])] + [sep.join([r[c] for c in columns]) for r in self._metrics])
        return "\n".join([sep.join([x for x in columns])] + [sep.join([self._format_field(r[c], fmt[c], ctype[c]) for c in columns]) for r in self._metrics])


class PicardHistMetrics(PicardMetrics):
    """Generic class to store metrics section from Picard Histogram
    Metrics reports. 

    In addition to metrics data the class also stores histogram
    values.

    Args:
      name (str): unique identifier for object, e.g. sample or unit name

      filename (str): file from which to collect metrics
      
      hist (list): list of lists, where the first entry holds the
        names of the values

      *args (list): if provided, must be a list of lists that upon
         initialization is converted to an OrderedDict. The first list
         is treated as a header.

    Returns:
      An instance of class PicardHistMetrics

    """

    def __init__(self, *args, identifier=None, filename=None, hist=None):
        # NB: __init__ should call super, but with current
        # implementation would require reading the metrics file twice!
        if filename is None and hist is None:
            raise ValueError("please provide argument hist when not reading from file")
        if filename is None and not args:
            raise ValueError("please supply either filename or args to instantiate class")
        self._set_vars(identifier, filename)
        if not self.filename is None and not args:
            (args, hist) = _read_picard_metrics(self.filename)
        self._metrics = DataFrame(*args)
        self._hist = DataFrame(*hist)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self.metrics._data]
        h = [self.hist.colnames] + [[row[c] for c in self.hist.colnames] for row in self.hist._data]
        return self.__class__(*a, identifier=self.id, filename=self.filename, hist=h)

    @property
    def hist(self):
        return self._hist


class AlignMetrics(PicardMetrics):
    _format = collections.OrderedDict([('CATEGORY', ('s', str)), ('TOTAL_READS', ('3.2h', int)), 
                                       ('PF_READS', ('3.2h', int)), ('PCT_PF_READS', ('3.2%', float)), 
                                       ('PF_NOISE_READS', ('3.2h', int)), ('PF_READS_ALIGNED', ('3.2h', int)), 
                                       ('PCT_PF_READS_ALIGNED', ('3.2%', float)), ('PF_ALIGNED_BASES', ('3.2h', int)), 
                                       ('PF_HQ_ALIGNED_READS', ('3.2h', int)), ('PF_HQ_ALIGNED_BASES', ('3.2h', int)), 
                                       ('PF_HQ_ALIGNED_Q20_BASES', ('3.2h', int)), ('PF_HQ_MEDIAN_MISMATCHES', ('', int)), 
                                       ('PF_MISMATCH_RATE', ('3.2f', float)), ('PF_HQ_ERROR_RATE', ('3.2f', float)), ('PF_INDEL_RATE', ('3.2f', float)), 
                                       ('MEAN_READ_LENGTH', ('3.2f', float)), ('READS_ALIGNED_IN_PAIRS', ('3.2h', int)), 
                                       ('PCT_READS_ALIGNED_IN_PAIRS', ('3.2%', float)), ('BAD_CYCLES', ('3.2h', int)), ('STRAND_BALANCE', ('3.2f', float)), 
                                       ('PCT_CHIMERAS', ('3.2%', float)), ('PCT_ADAPTER', ('3.2%', float)), ('SAMPLE', ('s', str)), 
                                       ('LIBRARY', ('s', str)), ('READ_GROUP', ('s', str))])

    def __init__(self, *args, identifier=None, filename=None):
        super(AlignMetrics, self).__init__(*args, identifier=identifier, filename=filename)

    def category(self, category="PAIR"):
        """Retrieve subset object with only one alignment category"""
        a = [self.metrics.colnames] + [[row[c] for c in self.metrics.colnames] for row in self._metrics if row['CATEGORY'] == category]
        return AlignMetrics(*a, filename=self.filename, identifier=self.id)

class InsertMetrics(PicardHistMetrics):
    _format = collections.OrderedDict([('MEDIAN_INSERT_SIZE', ('', int)), ('MEDIAN_ABSOLUTE_DEVIATION', ('', int)), 
                                       ('MIN_INSERT_SIZE', ('', int)), ('MAX_INSERT_SIZE', ('', int)), 
                                       ('MEAN_INSERT_SIZE', ('3.3f', float)), ('STANDARD_DEVIATION', ('3.3f', float)), 
                                       ('READ_PAIRS', ('3.2h', int)), ('PAIR_ORIENTATION', ('s', str)), 
                                       ('WIDTH_OF_10_PERCENT', ('', int)), ('WIDTH_OF_20_PERCENT', ('', int)),
                                       ('WIDTH_OF_30_PERCENT', ('', int)), ('WIDTH_OF_40_PERCENT', ('', int)), 
                                       ('WIDTH_OF_50_PERCENT', ('', int)), ('WIDTH_OF_60_PERCENT', ('', int)),
                                       ('WIDTH_OF_70_PERCENT', ('', int)), ('WIDTH_OF_80_PERCENT', ('', int)), 
                                       ('WIDTH_OF_90_PERCENT', ('', int)), ('WIDTH_OF_99_PERCENT', ('', int)),
                                       ('SAMPLE', ('s', str)), ('LIBRARY', ('s', str)), ('READ_GROUP', ('s', str))])
    def __init__(self, *args, identifier=None, filename=None, hist=None):
        super(InsertMetrics, self).__init__(*args, identifier=identifier, filename=filename, hist=hist)


class HsMetrics(PicardMetrics):
    _format = collections.OrderedDict([('BAIT_SET', ('s', str)), ('GENOME_SIZE', ('3.2h', int)), 
                                       ('BAIT_TERRITORY', ('3.2h', int)), ('TARGET_TERRITORY', ('3.2h', int)), 
                                       ('BAIT_DESIGN_EFFICIENCY', ('3.2f', float)), ('TOTAL_READS', ('3.2h', int)),
                                       ('PF_READS', ('3.2h', int)), ('PF_UNIQUE_READS', ('3.2h', int)), ('PCT_PF_READS', ('3.2%', float)), 
                                       ('PCT_PF_UQ_READS', ('3.2%', float)), ('PF_UQ_READS_ALIGNED', ('3.2h', int)), 
                                       ('PCT_PF_UQ_READS_ALIGNED', ('3.2%', float)), ('PF_UQ_BASES_ALIGNED', ('3.2h', int)), 
                                       ('ON_BAIT_BASES', ('3.2h', int)), ('NEAR_BAIT_BASES', ('3.2h', int)), ('OFF_BAIT_BASES', ('3.2h', int)),
                                       ('ON_TARGET_BASES', ('3.2h', int)), ('PCT_SELECTED_BASES', ('3.2%', float)), ('PCT_OFF_BAIT', ('3.2%', float)),
                                       ('ON_BAIT_VS_SELECTED', ('3.2f', float)), ('MEAN_BAIT_COVERAGE', ('3.2f', float)), 
                                       ('MEAN_TARGET_COVERAGE', ('3.2f', float)), ('PCT_USABLE_BASES_ON_BAIT', ('3.2%', float)),
                                       ('PCT_USABLE_BASES_ON_TARGET', ('3.2%', float)), ('FOLD_ENRICHMENT', ('3.2f', float)), 
                                       ('ZERO_CVG_TARGETS_PCT', ('3.2%', float)), ('FOLD_80_BASE_PENALTY', ('3.2f', float)), 
                                       ('PCT_TARGET_BASES_2X', ('3.2%', float)), ('PCT_TARGET_BASES_10X', ('3.2%', float)),
                                       ('PCT_TARGET_BASES_20X', ('3.2%', float)), ('PCT_TARGET_BASES_30X', ('3.2%', float)), 
                                       ('PCT_TARGET_BASES_40X', ('3.2%', float)), ('PCT_TARGET_BASES_50X', ('3.2%', float)), 
                                       ('PCT_TARGET_BASES_100X', ('3.2%', float)), ('HS_LIBRARY_SIZE', ('3.2h', int)), ('HS_PENALTY_10X', ('3.2f', float)),
                                       ('HS_PENALTY_20X', ('3.2f', float)), ('HS_PENALTY_30X', ('3.2f', float)), ('HS_PENALTY_40X', ('3.2f', float)),
                                       ('HS_PENALTY_50X', ('3.2f', float)), ('HS_PENALTY_100X', ('3.2f', float)), ('AT_DROPOUT', ('3.2f', float)), 
                                       ('GC_DROPOUT', ('3.2f', float)), ('SAMPLE', ('s', str)), ('LIBRARY',  ('s', str)), ('READ_GROUP',  ('s', str))])

    def __init__(self, *args, identifier=None, filename=None):
        super(HsMetrics, self).__init__(*args, identifier=identifier, filename=filename)


class DuplicationMetrics(PicardHistMetrics):
    _format = collections.OrderedDict([('LIBRARY', ('s', str)), ('UNPAIRED_READS_EXAMINED', ('3.2h', int)), 
                                       ('READ_PAIRS_EXAMINED', ('3.2h', int)), ('UNMAPPED_READS', ('3.2h', int)),
                                       ('UNPAIRED_READ_DUPLICATES', ('3.2h', int)), ('READ_PAIR_DUPLICATES', ('3.2h', int)), 
                                       ('READ_PAIR_OPTICAL_DUPLICATES', ('3.2f', float)), 
                                       ('PERCENT_DUPLICATION', ('3.2%', float)), ('ESTIMATED_LIBRARY_SIZE', ('3.2h', int))])

    def __init__(self, *args, identifier=None, filename=None, hist=None):
        super(DuplicationMetrics, self).__init__(*args, identifier=identifier, filename=filename, hist=hist)


def combine_metrics(metrics, mergename = "PicardMetrics_merge", uniquify=False):
    """Convert list of metrics objects to PicardMetrics object

    Args:
      metrics (list of tuples): list of metrics objects
      mergename (str): id to give to merged objet
      uniquify (bool): convert columns to unique names, appending .# where # is the count of the duplicated column

    Returns:
      new PicardMetrics object with format specifications set for summary operations
    """
    nrows = set([nr for sublist in  [[m.metrics.dim[0] for m in mtup] for mtup in metrics] for nr in sublist])
    if len(set(nrows)) > 1:
        raise ValueError("not all metrics of equal length; most probably you need to subset an AlignMetrics class to one category")
    ncols = set([nr for sublist in  [[m.metrics.dim[1] for m in mtup] for mtup in metrics] for nr in sublist])
    if len(set(ncols)) > len(metrics[0]):
        raise ValueError("not all metrics tuples have same set of columns; refusing to merge")
    colnames = [c for sublist in [m.metrics.colnames for m in metrics[0]] for c in sublist]
    args = [_make_unique(colnames)] if uniquify else [colnames]
    for mtup in metrics:
        rowargs = []
        fmtlist = []
        for m in mtup:
            rowargs += [m.metrics.data[0][k] for k in m.metrics.colnames]
            fmtlist += [(k,v) for k,v in m._format.items()]
        args.append(rowargs)
    pm = PicardMetrics(*args, identifier=mergename)
    pm._format = collections.OrderedDict(fmtlist)
    if uniquify:
        for c in args[0]:
            m = re.match("(.*)\.[0-9]+$", c)
            if m:
                pm._format[c] = pm._format[m.group(1)]
    return (pm)


