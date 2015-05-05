# Copyright (C) 2015 by Per Unneberg
import os
import re
import math
import pandas as pd
import numpy as np
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, gridplot
from bokeh.palettes import brewer
from snakemakelib.log import LoggerManager
from snakemakelib.bokeh.plot import lineplot, scatterplot2

smllogger = LoggerManager().getLogger(__name__)

# Picard has percentage columns that actually report fractions...
pct_re = re.compile(r'(PERCENT|PCT)')

def collect_picard_qc_results(inputfiles, samples):
    """Collect picard qc results"""
    first = True
    df_hist = df_met = None
    for (s, f) in zip(samples, inputfiles):
        smllogger.debug("parsing input file {f} for sample {s}".format(f=f, s=s))
        try:
            with open(f) as fh:
                data = [x.strip("\n").split("\t") for x in fh.readlines() if not x.strip() == ""]
            indices = list((i for i,val in enumerate(data) if val[0].startswith("## METRICS CLASS") or val[0].startswith("## HISTOGRAM")))
            if len(indices) == 1:
                indices += [len(data)]
            met_tmp = pd.DataFrame (data[indices[0]+2:indices[1]])
            met_tmp.columns = data[indices[0] + 1]
            pct_columns = [x for x in met_tmp.columns if pct_re.search(x)]
            for x in pct_columns:
                met_tmp[x] = met_tmp[x].astype(float) * 100.0
            met_tmp["Sample"] = s
            if indices[1] != len(data):
                hist_tmp = pd.DataFrame (data[indices[1]+2:])
                hist_tmp.columns = data[indices[1] + 1]
                hist_tmp["Sample"] = s
        except:
            smllogger.warn("pandas reading table {f} failed".format(f=f))
        try:
            if first:
                df_met = met_tmp
                if indices[1] != len(data):
                    df_hist = hist_tmp
                first = False
            else:
                df_met = df_met.append(met_tmp)
                if indices[1] != len(data):
                    df_hist = df_hist.append(hist_tmp)
        except:
                smllogger.warn("failed to append data")
    return (df_met, df_hist)


class Metrics(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super(Metrics, self).__init__(*args, **kwargs)
        self._metadata = {'type': 'metrics'}
        self.plots = []
        self._label = str(type(self)).split(".")[-1].replace("'>", "")

    @property
    def label(self):
        return self._label
        
    def plot_metrics(self, **kwargs):
        plist = []
        for (x,y,g) in self.plots:
            fig = scatterplot2(x=x, y=y, df=self, groups=g,
                               xaxis = {'axis_label' : x, 'major_label_orientation' : np.pi/3, 'axis_label_text_font_size' : '10pt'},
                              yaxis = {'axis_label' : y, 'major_label_orientation' : np.pi/3, 'axis_label_text_font_size' : '10pt'},
                               title="",
                              **kwargs)
            plist.append(fig)
        return plist
        
class HistMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(HistMetrics, self).__init__(*args, **kwargs)
        self._metadata = {'type' : 'histogram'}

    def plot_hist(self, **kwargs):
        fig = lineplot(self, groups=["Sample"])
        return fig
    
class AlignMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(AlignMetrics, self).__init__(*args, **kwargs)
        self.plots = [('Sample', 'PCT_PF_READS_ALIGNED', 'CATEGORY')]

class InsertMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(InsertMetrics, self).__init__(*args, **kwargs)
        self.plots = [('Sample','MEAN_INSERT_SIZE', [])]

class InsertHist(HistMetrics):
    def __init__(self, *args, **kwargs):
        super(InsertHist, self).__init__(*args, **kwargs)
        
class DuplicationMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(DuplicationMetrics, self).__init__(*args, **kwargs)
        self.plots = [('Sample', 'PERCENT_DUPLICATION', [])]

class DuplicationHist(HistMetrics):
    def __init__(self, *args, **kwargs):
        super(DuplicationHist, self).__init__(*args, **kwargs)
        
def _read_metrics(infile):
    if infile.endswith(".align_metrics.metrics.csv"):
        return AlignMetrics(pd.read_csv(infile))
    elif infile.endswith(".align_metrics.hist.csv"):
        return None
    elif infile.endswith(".insert_metrics.metrics.csv"):
        return InsertMetrics(pd.read_csv(infile))
    elif infile.endswith(".insert_metrics.hist.csv"):
        return InsertHist(pd.read_csv(infile))
    elif infile.endswith(".dup_metrics.metrics.csv"):
        return DuplicationMetrics(pd.read_csv(infile))
    elif infile.endswith(".dup_metrics.hist.csv"):
        return DuplicationHist(pd.read_csv(infile))
    else:
        smllogger.warn("Unknown metrics type {f}; skipping".format(f=infile))
        return None
    
def make_picard_summary_plots(inputfiles):
    d = {}
    TOOLS="pan,box_zoom,wheel_zoom,box_select,lasso_select,reset,save,hover"
    for (metrics_file, hist_file) in zip(inputfiles[0::2], inputfiles[1::2]):
        df_met = _read_metrics(metrics_file)
        df_hist = _read_metrics(hist_file)
        p1 = df_met.plot_metrics(tools=TOOLS)
        if not df_hist is None:
            p2 = [df_hist.plot_hist(tools=TOOLS)]
        else:
            p2 = []
        key = os.path.splitext(metrics_file)[0]
        if not df_met.label in d:
            d[df_met.label] = {}
        d[df_met.label][key] = gridplot([p1 + p2])
    return d

        
            
