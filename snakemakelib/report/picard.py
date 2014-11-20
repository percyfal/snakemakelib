# Copyright (c) 2014 Per Unneberg
import os
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
    _fieldnames = []
    _format = collections.OrderedDict()

    def __init__(self, *args):
        reader = csv.DictReader([",".join([str(y) for y in x]) for x in args])
        self._fieldnames = reader.fieldnames
        self._metrics = [collections.OrderedDict([(k, row[k]) for k in self._fieldnames]) for row in reader]

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
        return __class__(*a)

    @property
    def fieldnames(self):
        return self._fieldnames
    
    def summary(self, columns=None, fmt = None, ctype = None, sep="\t"):
        columns = list(self._format.keys()) if columns is None else columns
        fmt = {k:v[0] for (k,v) in list(self._format.items()) if k in columns} if fmt is None else {k:v for (k,v) in zip(columns, fmt)}
        ctype = {k:v[1] for (k,v) in list(self._format.items()) if k in columns} if ctype is None else {k:v for (k,v) in zip(columns, ctype)}
        return "\n".join([sep.join([x for x in columns])] + [sep.join(["{{{}}}".format(fmt[c]).format(ctype[c](r[c])) for c in columns]) for r in self._metrics])

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
                                       ('LIBRARY', (':s', str)), ('READ_GROUP', (':s', str))])
    def __init__(self, *args):
        super(AlignMetrics, self).__init__(*args)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return __class__(*a)


class InsertMetrics(PicardMetrics):
    _format = collections.OrderedDict([('MEDIAN_INSERT_SIZE', ('', int)), ('MEDIAN_ABSOLUTE_DEVIATION', ('', int)), 
                                       ('MIN_INSERT_SIZE', ('', int)), ('MAX_INSERT_SIZE', ('', int)), 
                                       ('MEAN_INSERT_SIZE', (':3.3f', float)), ('STANDARD_DEVIATION', (':3.3f', float)), 
                                       ('READ_PAIRS', (':3.2E', int)), ('PAIR_ORIENTATION', (':s', str)), 
                                       ('WIDTH_OF_10_PERCENT', ('', int)), ('WIDTH_OF_20_PERCENT', ('', int)),
                                       ('WIDTH_OF_30_PERCENT', ('', int)), ('WIDTH_OF_40_PERCENT', ('', int)), 
                                       ('WIDTH_OF_50_PERCENT', ('', int)), ('WIDTH_OF_60_PERCENT', ('', int)),
                                       ('WIDTH_OF_70_PERCENT', ('', int)), ('WIDTH_OF_80_PERCENT', ('', int)), 
                                       ('WIDTH_OF_90_PERCENT', ('', int)), ('WIDTH_OF_99_PERCENT', ('', int)),
                                       ('SAMPLE', (':s', str)), ('LIBRARY', (':s', str)), ('READ_GROUP', (':s', str))])
    def __init__(self, *args):
        super(InsertMetrics, self).__init__(*args)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return __class__(*a)


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
                                       ('GC_DROPOUT', (':3.2f', float)), ('SAMPLE', (':s', str)), ('LIBRARY',  (':s', str)), ('READ_GROUP',  (':s', str))])

    def __init__(self, *args):
        super(HsMetrics, self).__init__(*args)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return __class__(*a)


class DuplicationMetrics(PicardMetrics):
    _format = collections.OrderedDict([('LIBRARY', (':s', str)), ('UNPAIRED_READS_EXAMINED', (':3.2E', int)), 
                                       ('READ_PAIRS_EXAMINED', (':3.2E', int)), ('UNMAPPED_READS', (':3.2E', int)),
                                       ('UNPAIRED_READ_DUPLICATES', (':3.2E', int)), ('READ_PAIR_DUPLICATES', (':3.2E', int)), 
                                       ('READ_PAIR_OPTICAL_DUPLICATES', (':3.2f', float)), 
                                       ('PERCENT_DUPLICATION', (':3.2f', float)), ('ESTIMATED_LIBRARY_SIZE', (':3.2E', int))])

    def __init__(self, *args):
        super(DuplicationMetrics, self).__init__(*args)

    def __getitem__(self, columns):
        a = [columns] + [[row[c] for c in columns] for row in self._metrics]
        return __class__(*a)

class PicardHist(object):
    _fieldnames = []
    def __init__(self, *args):
        self._fieldnames = args[0]
        self._hist = collections.OrderedDict([(y[0], y[1:]) for y in [[x[i] for x in args] for i in range(0, len(self._fieldnames))]])

    @property
    def hist(self):
        return self._hist

    def __str__(self):
        return str(self._hist)

    @property
    def fieldnames(self):
        return self._fieldnames

class PicardMetricsCollection(object):
    """PicardMetrics: class for reading/storing metrics from one picard
    metrics file. The collection refers to the fact that often the
    files contain more than one type of metrics.

    Args:
      pmid (str): unique identifier
      file (str): metrics file to read

    """
    def __init__(self, pmid, file):
        self._metrics = None
        self._hist = None
        self._file = file
        self._pmid = pmid
        (_, self._metrics_type) = (os.path.splitext(file))
        self._read_metrics()

    def _read_metrics(self):
        """Read metrics"""
        (self._metrics, self._hist) = EXTENSIONS[self._metrics_type][2](self._file)
        # Possibly make separate containers for these
        self._update_data()

    def _update_data(self):
        reader = csv.DictReader([",".join([str(y) for y in x]) for x in self._metrics])
        self._metrics = [collections.OrderedDict([(k, row[k]) for k in reader.fieldnames]) for row in reader]
        # FIXME: this could probably be done more efficiently with some list comprehension
        if not self._hist is None:
            ncols = len(self._hist[0])
            nrows = len(self._hist)
            d = {x:[] for x in self._hist[0]}
            i = 0
            for c in self._hist[0]:
                d[c] = [x[i] for x in self._hist[1:]]
                i += 1
            self._hist = d

    def metrics(self, columns = None):
        if not columns is None:
            if not set(columns).issubset(set(self._metrics[0].keys())):
                raise KeyError("No such keys {} in {}".format(set(columns).difference(set(self._metrics[0].keys())), __name__))
            return [{k:v for k,v in x.items() if k in columns} for x in self._metrics]
        return self._metrics

    def hist(self):
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

    def __str__(self):
        return str(self._metrics)

    def __add__(self, other):
        pass


class AlignMetricsCollection(PicardMetricsCollection):
    def __init__(self, pmid, file):
        super(AlignMetricsCollection, self).__init__(pmid, file)

    def metrics(self, category = ['FIRST_OF_PAIR', 'SECOND_OF_PAIR', 'PAIR'], columns = None):
        metrics = [x for x in self._metrics if x['CATEGORY'] in category]
        if columns:
            if not set(columns).issubset(set(metrics[0].keys())):
                raise KeyError("No such keys {} in {}".format(set(columns).difference(set(metrics[0].keys())), __name__))
            return [{k:v for k,v in x.items() if k in columns and x['CATEGORY'] in category} for x in metrics]
        #return [{k:v for k,v in x.items()} for x in metrics]
        return metrics

class InsertMetricsCollection(PicardMetricsCollection):
    def __init__(self, pmid, file):
        super(InsertMetricsCollection, self).__init__(pmid, file)

class DuplicationMetricsCollection(PicardMetricsCollection):
    def __init__(self, pmid, file):
        super(DuplicationMetricsCollection, self).__init__(pmid, file)
        
class HsMetricsCollection(PicardMetricsCollection):
    def __init__(self, pmid, file):
        super(HsMetricsCollection, self).__init__(pmid, file)

class PicardMetricsSummary(object):
    """Combine different picard metrics to one summary.

    """

    def __init__(self, alnmetrics=None, hsmetrics=None, insertmetrics=None, dupmetrics=None):
        self.alnmetrics=alnmetrics
        self.hsmetrics=hsmetrics
        self.insertmetrics=insertmetrics
        self.dupmetrics=dupmetrics

    # def __str__(self):
    #     pass

