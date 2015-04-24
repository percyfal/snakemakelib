# Copyright (C) 2015 by Per Unneberg
import os
import re
import io
import pandas as pd
import jinja2
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, BoxSelectTool
from bokeh.models.widgets import VBox, HBox, TableColumn, DataTable
from bokeh.plotting import figure, output_file, show, gridplot
from bokeh.palettes import brewer
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

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

def make_picard_summary_plots(inputfiles):
    for (metrics_file, hist_file) in zip(inputfiles[0::2], inputfiles[1::2]):
        print (metrics_file, hist_file)
        df_met = pd.read_csv(metrics_file)
        df_hist = pd.read_csv(hist_file)
        if df_hist.endswith("has no histogram data"):
            df_hist = None

