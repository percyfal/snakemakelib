# Copyright (C) 2015 by Per Unneberg
import os
import re
import io
import pandas as pd
import jinja2
import numpy as np
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


                
