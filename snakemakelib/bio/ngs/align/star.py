# Copyright (C) 2015 by Per Unneberg
import pandas as pd
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, BoxSelectTool
from bokeh.models.widgets import VBox, HBox, TableColumn, DataTable
from bokeh.plotting import gridplot
from bokeh.palettes import brewer
from snakemake.report import data_uri
from snakemakelib.results import Results
from snakemakelib.report.utils import recast, trim_header
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)


class Star(Results):
    """Star: container class for star results"""
    _keys = ['align']

    def __init__(self, *args, **kw):
        self['align'] = None
        super(Star, self).__init__(*args, **kw)

    def _collect_results(self):
        smllogger.info("collecting results")
        df = None
        for (f, s) in zip(self._inputfiles, self._samples):
            smllogger.debug("Reading input file {f} for sample {s}".format(f=f, s=s))
            df_tmp = pd.read_table(f, sep="|",
                                   names=["name", "value"],
                                   engine="python", skiprows=[7, 22, 27])
            d = {trim_header(x, underscore=True, percent=True): recast(y)
                 for (x, y) in zip(df_tmp["name"], df_tmp["value"])}
            if df is None:
                df = pd.DataFrame(data=d, index=pd.Index([s], name="Sample"))
            else:
                df = df.append(pd.DataFrame(data=d, index=pd.Index([s], name="Sample")))
        df['mismatch_sum'] = df['Mismatch_rate_per_base__PCT'] +\
            df['Deletion_rate_per_base'] + df['Insertion_rate_per_base']
        df['PCT_of_reads_unmapped'] = df['PCT_of_reads_unmapped:_other'] +\
            df['PCT_of_reads_unmapped:_too_many_mismatches'] +\
            df['PCT_of_reads_unmapped:_too_short']
        self['align'] = df

