# Copyright (C) 2015 by Per Unneberg
import os
import re
import io
import pandas as pd
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import gridplot
from bokeh.palettes import brewer
from snakemake.report import data_uri
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)


def collect_rseqc_results(inputfiles, samples):
    """Collect rseqc results"""
    first = True
    df_rd = df_gc = None
    for (f, s) in zip(inputfiles, samples):
        smllogger.debug("parsing input file {f} for sample {s}".format(f=f, s=s))
        try:
            infile = os.path.join(os.path.dirname(f), "read_distribution.txt")
            with open(infile, "r") as fh:
                data = fh.readlines()
            d1 = io.StringIO("".join([re.sub(r"(\w) (\w)", r"\1_\2", x)
                                      for x in data[3:6]]))
            tmp = pd.read_table(d1, engine="python", sep="[ ]+", header=None)
            tmp.columns = ["Group", "Tag_count"]
            tmp["Total_bases"] = None
            tmp["Tags/Kb"] = None
            d2 = io.StringIO("".join(data[7:17]))
            df_rd_tmp = pd.read_table(d2, engine="python", sep="[ ]+")
            df_rd_tmp = df_rd_tmp.append(tmp)
            df_rd_tmp["Sample"] = s
        except:
            smllogger.warn("pandas reading table {f} failed".format(f=infile))
        try:
            infile = os.path.join(os.path.dirname(f),
                                  "geneBody_coverage.geneBodyCoverage.txt")
            df_gc_tmp = pd.read_table(infile, engine="python", sep="\t")
            df_gc_tmp["Sample"] = s
        except:
            smllogger.warn("pandas reading table {f} failed".format(f=infile))
        try:
            if first:
                df_rd = df_rd_tmp
                df_gc = df_gc_tmp
                first = False
            else:
                df_rd = df_rd.append(df_rd_tmp)
                df_gc = df_gc.append(df_gc_tmp)
        except:
            smllogger.warn("failed to append data")
    df_rd = df_rd.set_index(["Group", "Sample"])
    df_gc = df_gc.set_index("Sample")
    df_rd = df_rd.astype(float)
    ExonMap = 100 * (df_rd.loc['CDS_Exons'] + df_rd.loc["3'UTR_Exons"] + df_rd.loc["5'UTR_Exons"]) / df_rd.loc['Total_Assigned_Tags']
    ExonMap["Group"] = "ExonMap"
    ExonMap = ExonMap.reset_index().set_index(["Group", "Sample"])
    df_rd = df_rd.append(ExonMap)
    df_gc["three_prime_map"] = 100.0 * df_gc.loc[:, "91":"100"].sum(axis=1) /\
                               df_gc.loc[:, "1":"100"].sum(axis=1)
    return {'rd': df_rd.sortlevel(1), 'gc': df_gc}


