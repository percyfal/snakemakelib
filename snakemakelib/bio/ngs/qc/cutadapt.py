# Copyright (C) 2015 by Per Unneberg
import re
import pandas as pd
import numpy as np
from bokeh.models import HoverTool
from snakemake.report import data_uri
from snakemakelib.log import LoggerManager
from snakemakelib.bokeh.plot import make_dotplot

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
    return df_ret


def make_cutadapt_summary_plot(inputfile):
    df_summary = pd.read_csv(inputfile)
    df_summary["read1_pct"] = 100.0 * df_summary["Read 1 with adapter"] /\
        df_summary["Total read pairs processed"]
    df_summary["read2_pct"] = 100.0 * df_summary["Read 2 with adapter"] /\
        df_summary["Total read pairs processed"]

    TOOLS = "pan,wheel_zoom,box_zoom,box_select,reset,save,hover"
    fig = make_dotplot(y=["read1_pct", "read2_pct"], df=df_summary,
                       groups=["sample", "run"],
                       tooltips=[{'type': HoverTool,
                                  'tips': [('Sample', '@sample'), ]}],
                       plot_width=400, plot_height=400, tools=TOOLS,
                       title="Cutadapt metrics",
                       title_text_font_size='12pt',
                       xaxis={'axis_label': "sample",
                              'major_label_orientation': np.pi/3,
                              'axis_label_text_font_size': '10pt'},
                       yaxis={'axis_label': "percent reads",
                              'major_label_orientation': 1,
                              'axis_label_text_font_size': '10pt'},
                       x_axis_type=None, y_axis_type="linear",
                       circle={'size': 10, 'alpha': 0.5,
                               'line_color': 'black'},
                       sidelegend=False,
                       y_range=[0, 105])

    return {'fig': fig, 'uri': data_uri(inputfile), 'file': inputfile}
