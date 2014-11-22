# Copyright (c) 2014 Per Unneberg
import os
import sys
import re
import csv
import collections

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

# For now: extension maps to tuple (label, description). Label should
# be reused for analysis definitions
EXTENSIONS={'.align_metrics':('align', 'alignment', _read_picard_metrics),
            '.hs_metrics':('hs', 'hybrid selection', _read_picard_metrics),
            '.dup_metrics':('dup', 'duplication metrics', _read_picard_metrics),
            '.insert_metrics':('insert', 'insert size', _read_picard_metrics),
            '.eval_metrics':('eval', 'snp evaluation', _raw)
            }

class PicardMetrics(object):
    """Generic class to store metrics section from Picard Metrics reports.
    See also class PicardHistMetrics for reports that provide metrics
    and histogram information.


    Args:
      name (str): unique identifier for object, e.g. sample or unit name
      filename (str): file from which to collect metrics
      *args: if provided, must be a list of lists that upon initialization is converted to an OrderedDict. The first list is treated as a header.

    Returns:
      An instance of class PicardMetrics
    """
    _format = collections.OrderedDict()

    def __init__(self, *args, filename=None, identifier=None):
        if filename is None and not args:
            raise ValueError("please supply either filename or args to instantiate class")
        self._set_vars(identifier, filename)
        if not self.filename is None and not args:
            (args, _) = _read_picard_metrics(self.filename)
        self._set_metrics(args)

    def _set_vars(self, identifier, filename):
        self._id = identifier if not identifier is None else filename
        self._filename = filename

    def _set_metrics(self, args):
        reader = csv.DictReader([",".join([str(y) for y in x]) for x in args])
        self._metrics = [collections.OrderedDict([(k, row[k]) for k in reader.fieldnames] + [("ID", self.id)]) for row in reader]
        self._fieldnames = reader.fieldnames + ["ID"]

    def __str__(self):
        return str(self._metrics)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if len(self._metrics) > self.index:
            self.index += 1
            return self._metrics[self.index - 1]
        else:
            raise StopIteration

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return PicardMetrics(*a, filename=self.filename, identifier=self.id)

    def __add__(self, other):
        if not self.fieldnames == other.fieldnames:
            raise TypeError("fieldnames differ between {id1} and {id2}: {fn1} != {fn2}; cannot merge objects with different fieldnames".format(id1=self.id, id2=other.id, fn1=self.fieldnames, fn2=other.fieldnames))
        a = [self.fieldnames] + [[row[c] for c in self.fieldnames] for row in self._metrics] + [[row[c] for c in other.fieldnames] for row in other.metrics]
        #return PicardMetrics(*a, filename=",".join([self.filename, other.filename]), identifier=",".join([self.id, other.id]))
        return PicardMetrics(*a, filename=",".join([self.filename, other.filename]))#, identifier=",".join([self.id, other.id]))

    @property
    def fieldnames(self):
        return self._fieldnames

    @property
    def metrics(self):
        return self._metrics
    
    @property 
    def id(self):
        return self._id

    @property
    def filename(self):
        return self._filename

    def summary(self, fmt = None, ctype = None, sep="\t"):
        columns = self.fieldnames
        fmt = {k:v[0] for (k,v) in list(self._format.items()) if k in columns} if fmt is None else {k:v for (k,v) in zip(columns, fmt)}
        ctype = {k:v[1] for (k,v) in list(self._format.items()) if k in columns} if ctype is None else {k:v for (k,v) in zip(columns, ctype)}
        if not fmt:
            raise ValueError("No format defined for {cls}; please used derived subclass".format(cls=__class__))
        return "\n".join([sep.join([x for x in columns])] + [sep.join(["{{{}}}".format(fmt[c]).format(ctype[c](r[c])) for c in columns]) for r in self._metrics])


class PicardHistMetrics(PicardMetrics):
    """Generic class to store metrics section from Picard Histogram
    Metrics reports. 

    In addition to metrics data the class also stores histogram
    values. Nevertheless, iterations and slices of the class will only
    be applied on the metrics tables as it makes little sense to slice
    or iterate over the histograms.

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
        self._set_metrics(args)
        self._histfieldnames = hist[0]
        self._hist = collections.OrderedDict([(y[0], y[1:]) for y in [[x[i] for x in hist] for i in range(0, len(self._histfieldnames))]])

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return PicardHistMetrics(*a, filename=self.filename, identifier=self.id, hist=list(self.hist))
        
    @property
    def hist(self):
        return self._hist
    
    @property
    def histfieldnames(self):
        return self._histfieldnames

class AlignMetrics(PicardMetrics):
    _format = collections.OrderedDict([('CATEGORY', (':s', str)), ('TOTAL_READS', (':3.2E', int)), 
                                       ('PF_READS', (':3.2E', int)), ('PCT_PF_READS', (':3.2f', float)), 
                                       ('PF_NOISE_READS', (':3.2E', int)), ('PF_READS_ALIGNED', (':3.2E', int)), 
                                       ('PCT_PF_READS_ALIGNED', (':3.2f', float)), ('PF_ALIGNED_BASES', (':3.2E', int)), 
                                       ('PF_HQ_ALIGNED_READS', (':3.2E', int)), ('PF_HQ_ALIGNED_BASES', (':3.2E', int)), 
                                       ('PF_HQ_ALIGNED_Q20_BASES', (':3.2E', int)), ('PF_HQ_MEDIAN_MISMATCHES', (':3.2E', int)), 
                                       ('PF_MISMATCH_RATE', (':3.2f', float)), ('PF_HQ_ERROR_RATE', (':3.2f', float)), ('PF_INDEL_RATE', (':3.2f', float)), 
                                       ('MEAN_READ_LENGTH', (':3.2f', float)), ('READS_ALIGNED_IN_PAIRS', (':3.2E', int)), 
                                       ('PCT_READS_ALIGNED_IN_PAIRS', (':3.2f', float)), ('BAD_CYCLES', (':3.2E', int)), ('STRAND_BALANCE', (':3.2f', float)), 
                                       ('PCT_CHIMERAS', (':3.2f', float)), ('PCT_ADAPTER', (':3.2f', float)), ('SAMPLE', (':s', str)), 
                                       ('LIBRARY', (':s', str)), ('READ_GROUP', (':s', str)), ('ID', (':s', str))])

    def __init__(self, *args, identifier=None, filename=None):
        super(AlignMetrics, self).__init__(*args, identifier=identifier, filename=filename)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return AlignMetrics(*a, filename=self.filename, identifier=self.id)

class InsertMetrics(PicardHistMetrics):
    _format = collections.OrderedDict([('MEDIAN_INSERT_SIZE', ('', int)), ('MEDIAN_ABSOLUTE_DEVIATION', ('', int)), 
                                       ('MIN_INSERT_SIZE', ('', int)), ('MAX_INSERT_SIZE', ('', int)), 
                                       ('MEAN_INSERT_SIZE', (':3.3f', float)), ('STANDARD_DEVIATION', (':3.3f', float)), 
                                       ('READ_PAIRS', (':3.2E', int)), ('PAIR_ORIENTATION', (':s', str)), 
                                       ('WIDTH_OF_10_PERCENT', ('', int)), ('WIDTH_OF_20_PERCENT', ('', int)),
                                       ('WIDTH_OF_30_PERCENT', ('', int)), ('WIDTH_OF_40_PERCENT', ('', int)), 
                                       ('WIDTH_OF_50_PERCENT', ('', int)), ('WIDTH_OF_60_PERCENT', ('', int)),
                                       ('WIDTH_OF_70_PERCENT', ('', int)), ('WIDTH_OF_80_PERCENT', ('', int)), 
                                       ('WIDTH_OF_90_PERCENT', ('', int)), ('WIDTH_OF_99_PERCENT', ('', int)),
                                       ('SAMPLE', (':s', str)), ('LIBRARY', (':s', str)), ('READ_GROUP', (':s', str)), ('ID', (':s', str))])
    def __init__(self, *args, identifier=None, filename=None, hist=None):
        super(InsertMetrics, self).__init__(*args, identifier=identifier, filename=filename, hist=hist)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return InsertMetrics(*a, filename=self.filename, identifier=self.id, hist=list(self.hist))


class HsMetrics(PicardMetrics):
    _format = collections.OrderedDict([('BAIT_SET', (':s', str)), ('GENOME_SIZE', (':3.2E', int)), 
                                       ('BAIT_TERRITORY', (':3.2E', int)), ('TARGET_TERRITORY', (':3.2E', int)), 
                                       ('BAIT_DESIGN_EFFICIENCY', (':3.2f', float)), ('TOTAL_READS', (':3.2E', int)),
                                       ('PF_READS', (':3.2E', int)), ('PF_UNIQUE_READS', (':3.2E', int)), ('PCT_PF_READS', (':3.2f', float)), 
                                       ('PCT_PF_UQ_READS', (':3.2f', float)), ('PF_UQ_READS_ALIGNED', (':3.2E', int)), 
                                       ('PCT_PF_UQ_READS_ALIGNED', (':3.2f', float)), ('PF_UQ_BASES_ALIGNED', (':3.2E', int)), 
                                       ('ON_BAIT_BASES', (':3.2E', int)), ('NEAR_BAIT_BASES', (':3.2E', int)), ('OFF_BAIT_BASES', (':3.2E', int)),
                                       ('ON_TARGET_BASES', (':3.2E', int)), ('PCT_SELECTED_BASES', (':3.2f', float)), ('PCT_OFF_BAIT', (':3.2f', float)),
                                       ('ON_BAIT_VS_SELECTED', (':3.2f', float)), ('MEAN_BAIT_COVERAGE', (':3.2f', float)), 
                                       ('MEAN_TARGET_COVERAGE', (':3.2f', float)), ('PCT_USABLE_BASES_ON_BAIT', (':3.2f', float)),
                                       ('PCT_USABLE_BASES_ON_TARGET', (':3.2f', float)), ('FOLD_ENRICHMENT', (':3.2f', float)), 
                                       ('ZERO_CVG_TARGETS_PCT', (':3.2f', float)), ('FOLD_80_BASE_PENALTY', (':3.2f', float)), 
                                       ('PCT_TARGET_BASES_2X', (':3.2f', float)), ('PCT_TARGET_BASES_10X', (':3.2f', float)),
                                       ('PCT_TARGET_BASES_20X', (':3.2f', float)), ('PCT_TARGET_BASES_30X', (':3.2f', float)), 
                                       ('PCT_TARGET_BASES_40X', (':3.2f', float)), ('PCT_TARGET_BASES_50X', (':3.2f', float)), 
                                       ('PCT_TARGET_BASES_100X', (':3.2f', float)), ('HS_LIBRARY_SIZE', (':3.2E', int)), ('HS_PENALTY_10X', (':3.2f', float)),
                                       ('HS_PENALTY_20X', (':3.2f', float)), ('HS_PENALTY_30X', (':3.2f', float)), ('HS_PENALTY_40X', (':3.2f', float)),
                                       ('HS_PENALTY_50X', (':3.2f', float)), ('HS_PENALTY_100X', (':3.2f', float)), ('AT_DROPOUT', (':3.2f', float)), 
                                       ('GC_DROPOUT', (':3.2f', float)), ('SAMPLE', (':s', str)), ('LIBRARY',  (':s', str)), ('READ_GROUP',  (':s', str)), ('ID', (':s', str))])

    def __init__(self, *args, identifier=None, filename=None):
        super(HsMetrics, self).__init__(*args, identifier=identifier, filename=filename)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return HsMetrics(*a)


class DuplicationMetrics(PicardHistMetrics):
    _format = collections.OrderedDict([('LIBRARY', (':s', str)), ('UNPAIRED_READS_EXAMINED', (':3.2E', int)), 
                                       ('READ_PAIRS_EXAMINED', (':3.2E', int)), ('UNMAPPED_READS', (':3.2E', int)),
                                       ('UNPAIRED_READ_DUPLICATES', (':3.2E', int)), ('READ_PAIR_DUPLICATES', (':3.2E', int)), 
                                       ('READ_PAIR_OPTICAL_DUPLICATES', (':3.2f', float)), 
                                       ('PERCENT_DUPLICATION', (':3.2f', float)), ('ESTIMATED_LIBRARY_SIZE', (':3.2E', int)), ('ID', (':s', str))])

    def __init__(self, *args, identifier=None, filename=None, hist=None):
        super(DuplicationMetrics, self).__init__(*args, identifier=identifier, filename=filename, hist=hist)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return DuplicationMetrics(*a, filename=self.filename, identifier=self.id, hist=list(self.hist))


class PicardMetricsSummary(object):
    """Combine different picard metrics to one summary.

    Takes AlignMetrics etc as input and outputs the *metrics*
    summaries only. If all alnmetrics lines chosen, replicate the
    others three times.

    """

    def __init__(self, alnmetrics=None, hsmetrics=None, insertmetrics=None, dupmetrics=None):
        self.alnmetrics=alnmetrics
        self.hsmetrics=hsmetrics
        self.insertmetrics=insertmetrics
        self.dupmetrics=dupmetrics

    # def __str__(self):
    #     pass

