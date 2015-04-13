# Copyright (C) 2015 by Per Unneberg
import pandas as pd
import jinja2
import numpy as np
from snakemakelib.report.utils import recast, trim_header

def collect_star_alignment_results(input, samples):
    """Collect star alignment results"""
    df = None
    for (f, s) in zip(input, samples):
        df_tmp = pd.read_table(f, sep="|", names=["name", "value"], engine="python", skiprows=[7,22,27])
        d = {trim_header(x, percent=True):recast(y) for (x,y) in zip(df_tmp["name"], df_tmp["value"])}
        if df is None:
            df = pd.DataFrame(d, index=[s])
        else:
            df = df.append(pd.DataFrame(d, index=[s]))
    return df
