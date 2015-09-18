# Copyright (C) 2015 by Per Unneberg
import re
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import HoverTool
from snakemake.report import data_uri
from snakemakelib.log import LoggerManager
from bokehutils.mgeom import mdotplot
from bokehutils.axes import xaxis, yaxis

smllogger = LoggerManager().getLogger(__name__)


def collect_cutadapt_qc_results(inputfiles, sampleruns):
    """Collect cutadapt metrics"""
    first = True
    df = None
    for (sr, f) in zip(sampleruns, inputfiles):
        (sample, run) = sr
        smllogger.debug("parsing input file {f} for sample {s}, run {r}".format(f=f, s=sample, r=run))
        try:
            with open(f) as fh:
                data = [x.strip("\n").split("\t") for x in fh.readlines()
                        if not x.strip() == ""]
            indices = list((i for i, val in enumerate(data)
                            if val[0].startswith("===")))
            summary = (list([
                re.sub(r'([ ]?\(\d+\.\d+%\)|,|[ ]?bp| $)', r'', z.strip(" "))
                for z in x
            ]
                            for x in list(y.split(":")
                                          for x in data[indices[0]+1:indices[1]]
                                          for y in x)))
            df_summary = pd.DataFrame(
                summary, columns=["statistic", "count"])
            df_summary["sample"] = sample
            df_summary["run"] = run
        except:
            smllogger.warn("pandas reading metrics file {f} failed".format(f=f))
        try:
            if first:
                df = df_summary
                first = False
            else:
                df = df.append(df_summary)
        except:
            smllogger.warn("failed to append data")
    df[['count']] = df[['count']].astype(int)
    df_ret = df.pivot_table(
        columns=["statistic"], values=["count"], index=["sample", "run"])
    df_ret.columns = df_ret.columns.droplevel()
    df_ret["read1_pct"] = 100.0 * df_ret["Read 1 with adapter"] /\
        df_ret["Total read pairs processed"]
    df_ret["read2_pct"] = 100.0 * df_ret["Read 2 with adapter"] /\
        df_ret["Total read pairs processed"]

    return df_ret

def make_cutadapt_summary_plot(inputfile):
    df_summary = pd.read_csv(inputfile)
    df_summary["sample"] = df_summary["sample"].astype("str")
    TOOLS = "pan,wheel_zoom,box_zoom,box_select,reset,save"
    fig = figure(tools=TOOLS, width=400, height=400,
                 x_range=list(set(df_summary["sample"])),
                 y_range=[0, 105], title="Cutadapt metrics",
                 title_text_font_size='12pt')
    mdotplot(fig, x="sample", y=["read1_pct", "read2_pct"],
             df=df_summary, size=10, alpha=0.5)
    xaxis(fig, axis_label="sample",
          major_label_orientation=np.pi/3,
          axis_label_text_font_size='10pt')
    yaxis(fig, axis_label="percent reads",
          major_label_orientation=1,
          axis_label_text_font_size='10pt')
    return {'fig': fig, 'uri': data_uri(inputfile), 'file': inputfile}
