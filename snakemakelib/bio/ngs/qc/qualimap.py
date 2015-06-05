# Copyright (C) 2015 by Per Unneberg
import os
import re
import io
import pandas as pd
import jinja2
import numpy as np
from bokeh.plotting import figure
from bokeh.models import GridPlot
from snakemakelib.utils import SmlTemplateEnv
from snakemake.report import data_uri
from snakemakelib.results import Results
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

class Qualimap(Results):
    _keys = ['coverage_per_contig']
    
    def __init__(self, *args, **kw):
        super(Qualimap, self).__init__(*args, **kw)

    def _collect_results(self):
        smllogger.info("Collecting results")
        first = True
        for (f, s) in zip(self._inputfiles, self._samples):
            data = self.load_lines(f)
            df_tmp = self.parse_data(data, rs=("Coverage per contig", None), skip=2, split=True, columns=["chr", "chrlen", "mapped_bases", "mean_coverage", "sd"], dtype=float)
            df_tmp["Sample"] = s
            try:
                if first:
                    df = df_tmp
                    first = False
                else:
                    df.append(df_tmp)
            except:
                smllogger.warn("failed to append data to coverage_per_contig dataframe")
        df['chrlen_percent'] = 100 * df['chrlen']/sum(df['chrlen'])
        df['mapped_bases_percent'] = 100 * df['mapped_bases']/sum(df['mapped_bases'])
        self['coverage_per_contig'] = df

def make_qualimap_plots(coverage_per_contig = None,
                        **kwargs):
    """Make qualimap summary plots"""
    df_all = pd.read_csv(coverage_per_contig, index_col=0)
    samples = list(df_all.Sample)
    plist = []
    for s in samples:
        df = df_all[df_all["Sample"] == s]
        p = figure(title="plot", 
                   x_range=[-1, max(df['chrlen_percent']) * 1.1], 
                   y_range = [-1, max(df['mapped_bases_percent']) * 1.1])
        p.title = "Sample {}".format(s)
        p.text(df['chrlen_percent'], df['mapped_bases_percent'], text=df['chr'])
        m = max(df['chrlen_percent'] + df['mapped_bases_percent'])
        p.line(x=[0,m], y=[0,m])
        plist.append(p)
    children = [plist[i:i+3] for i in range(0, len(plist), 3)]
    return {'fig' : GridPlot(children=children, title="Fraction mapped bases vs fraction chrlen"),
            'uri' : data_uri(coverage_per_contig),
            'file' : coverage_per_contig}
    

    
