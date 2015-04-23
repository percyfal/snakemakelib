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
            if f.endswith(".insert_metrics"):
                indices = list((i for i,val in enumerate(data) if val[0].startswith("## METRICS CLASS") or val[0].startswith("## HISTOGRAM")))
                met_tmp = pd.DataFrame (data[indices[0]+2:indices[1]])
                met_tmp.columns = data[indices[0] + 1]
                print ("sample: ", s)
                met_tmp.index = s
                #print (s, met_tmp)
                print (s)
                
                hist_tmp = pd.DataFrame (data[indices[1]+2:], index=s)
                hist_tmp.columns = data[indices[1] + 1]
        except:
            smllogger.warn("pandas reading table {f} failed".format(f=f))
        try:
            if first:
                df_met = met_tmp
                df_hist = hist_tmp
                first = False
            else:
                df_met = df_met.append(met_tmp)
                df_hist = df_hist.append(hist_temp)
        except:
                smllogger.warn("failed to append data")
    return {'metrics' : df_met, 'hist' : df_hist}

def make_picard_summary_plots(d):
    pass
# def make_picard_summary_plots(df_rd, df_gc, do_qc=True, min_exonmap=60.0, max_three_prime_map=10.0):
#     """Make picard summary plots"""
#     pass
