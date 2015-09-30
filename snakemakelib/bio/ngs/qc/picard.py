# Copyright (C) 2015 by Per Unneberg
import os
import re
import pandas as pd
import numpy as np
from bokeh.charts import Scatter
from bokeh.plotting import figure, gridplot
from bokehutils.geom import dotplot, lines
from bokehutils.mgeom import mlines
from bokehutils.facet import facet_grid
from bokehutils.axes import xaxis, yaxis
from bokehutils.color import colorbrewer
from snakemake.report import data_uri
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

# Picard has percentage columns that actually report fractions...
pct_re = re.compile(r'(PERCENT|PCT)')

def collect_picard_qc_results(inputfiles, samples):
    """Collect picard qc results"""
    first = True
    df_hist = df_met = None
    for (s, f) in zip(samples, inputfiles):
        smllogger.debug("parsing input file {f} for sample {s}".format(
            f=f, s=s))
        try:
            with open(f) as fh:
                data = [x.strip("\n").split("\t") for x in fh.readlines()
                        if not x.strip() == ""]
            indices = list((i for i, val in enumerate(data)
                            if val[0].startswith("## METRICS CLASS")
                            or val[0].startswith("## HISTOGRAM")))
            if len(indices) == 1:
                indices += [len(data)]
            met_tmp = pd.DataFrame(data[indices[0]+2:indices[1]])
            met_tmp.columns = data[indices[0] + 1]
            pct_columns = [x for x in met_tmp.columns if pct_re.search(x)]
            for x in pct_columns:
                met_tmp[x] = met_tmp[x].astype(float) * 100.0
            met_tmp["Sample"] = s
            if indices[1] != len(data):
                hist_tmp = pd.DataFrame(data[indices[1]+2:])
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
        self["Sample"] = self["Sample"].astype('str')
        self._metadata = {'type': 'metrics'}
        self.plots = []
        self._label = str(type(self)).split(".")[-1].replace("'>", "")
        self.ncol = 4
        self._stack = None
        self._kwargs = None

    @property
    def label(self):
        return self._label

    def plot_metrics(self, **kwargs):
        """Plot metrics wrapper

        Returns:
          plist (list): list of bokeh plot objects
        """
        plist = []
        for kw in self.plots:
            kwargs.update(kw['figure'])
            if self._stack is None:
                fig = figure(**kwargs)
                dotplot(fig, df=self, **kw['renderer'])
                xaxis(fig, **kw['xaxis'])
                yaxis(fig, **kw['yaxis'])
            else:
                fig = Scatter(self._stack, **self._kwargs)
            plist.append(fig)

        return plist


class HistMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(HistMetrics, self).__init__(*args, **kwargs)
        self._metadata = {'type': 'histogram'}
        self.kw = {}
        self.plot_hist = self.lines_hist

    def facet_grid_hist(self, **kwargs):
        """facet_grid hist wrapper.

        Returns:
          plist (list): list of bokeh plot objects
        """
        kwargs.update(self.kw)
        f = figure(**kwargs['figure'])
        mlines(f, df=self, **kwargs['renderer'])
        gp = facet_grid(f, df=self, ncol=self.ncol,
                        **kwargs['facet'])
        plist = [x for sublist in gp.children for x in sublist]
        return plist

    def lines_hist(self, **kwargs):
        """lines hist wrapper.

        Plot lines in one plot, with legend.

        Returns:
          plist (list): list of bokeh plot objects
        """
        kwargs.update(self.kw)
        f = figure(**kwargs['figure'])
        lines(f, df=self, **kwargs['renderer'])
        return [f]



class AlignMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(AlignMetrics, self).__init__(*args, **kwargs)
        self.plots = [{
            'figure': {
                'plot_width': 400, 'plot_height': 400,
                'title': 'Percent PF_READS aligned per sample',
                'y_range': [0, 100],
                'title_text_font_size': "10pt",
            },
            'renderer': {'size': 10, 'alpha': 0.5,
                         'x': 'Sample',
                         'y': 'PCT_PF_READS_ALIGNED',
                         #'groups': 'CATEGORY',
            },
            'xaxis': {'axis_label': 'Sample',
                      'major_label_orientation': np.pi/3,
                      'axis_label_text_font_size': '10pt'},
            'yaxis': {'axis_label': 'Percentage reads',
                      'major_label_orientation': np.pi/3,
                      'axis_label_text_font_size': '10pt'}}]


class InsertMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(InsertMetrics, self).__init__(*args, **kwargs)
        self.columns = ["ind"] + list(self.columns[1:])
        self.plots = [{
            'figure': {
                'plot_width': 400, 'plot_height': 400,
                'title': 'Mean insert size',
                'title_text_font_size': "10pt",
                'x_range': list(self.Sample),
            },
            'renderer': {
                'size': 10, 'alpha': 0.3,
                'x': 'Sample', 'y': 'MEAN_INSERT_SIZE',
            },
            'xaxis': {'axis_label': 'Sample',
                      'major_label_orientation': np.pi/3,
                      'axis_label_text_font_size': '10pt'
            },
            'yaxis': {'axis_label': 'Mean insert size',
                      'major_label_orientation': np.pi/3,
                      'axis_label_text_font_size': '10pt'}}]
        self._stack = self[["Sample", "ind", "MEAN_INSERT_SIZE"]].pivot_table(index=["Sample", "ind"]).stack().reset_index()
        self._stack.columns = ["Sample", "ind", "type", "MEAN_INSERT_SIZE"]
        self._kwargs = {'x': 'Sample', 'y': 'MEAN_INSERT_SIZE', 'color': 'ind',
                        'legend': 'top_right', 'title': "Picard metrics, mean insert size",
                        'ylabel': "mean insert size"}


class InsertHist(HistMetrics):
    def __init__(self, *args, **kwargs):
        super(InsertHist, self).__init__(*args, **kwargs)
        self.kw = {
            'figure': {
                'plot_width': 400, 'plot_height': 400,
                'title': "Insert size distribution",
                'title_text_font_size': "10pt",
            },
            'facet': {
                'groups': ["Sample"],
                'width': 300, 'height': 300,
                'share_x_range': True,
                'x': 'insert_size',
                'y':  [x for x in list(self.columns) if x.startswith("All_Reads")],
                'title_text_font_size': "12pt",
            },
            'renderer': {
                'x': 'insert_size',
                'y': [x for x in list(self.columns) if x.startswith("All_Reads")],
                'legend' : [x for x in list(self.columns) if x.startswith("All_Reads")],
                'line_width' : 2,
            },
            'xaxis': {'axis_label': 'Insert size',
                      'axis_label_text_font_size': '10pt'},
            'yaxis': {'axis_label': 'Count',
                      'axis_label_text_font_size': '10pt'}}
        self.kw['renderer']['color'] =  colorbrewer(datalen = len(self.kw['renderer']['y']))
        self.plot_hist = self.facet_grid_hist


class DuplicationMetrics(Metrics):
    def __init__(self, *args, **kwargs):
        super(DuplicationMetrics, self).__init__(*args, **kwargs)
        self.plots = [{
            'figure': {
                'plot_width': 400, 'plot_height': 400,
                'title': 'Percent duplication per sample',
                'title_text_font_size': "10pt",
                'y_range': [0, 100],
                'x_range': list(self.Sample),
            },
            'renderer': {
                'size': 10, 'alpha': 0.3,
                'x': 'Sample', 'y': 'PERCENT_DUPLICATION',
            },
            'xaxis': {'axis_label': 'Sample',
                      'major_label_orientation': np.pi/3,
                      'axis_label_text_font_size': '10pt'},
            'yaxis': {'axis_label': 'Percent duplication',
                      'major_label_orientation': np.pi/3,
                      'axis_label_text_font_size': '10pt'}}]
        self._stack = self[["Sample", "PERCENT_DUPLICATION"]].pivot_table(index="Sample").stack().reset_index([0,1])
        self._stack.columns = ["Sample", "PERCENT_DUPLICATION", "percent"]
        self._kwargs = {'x': 'Sample', 'y': 'percent',
                        'title': "Picard metrics, percent duplication"}

class DuplicationHist(HistMetrics):
    # see http://sourceforge.net/p/picard/wiki/Main_Page/
    # q-what-is-meaning-of-the-histogram-produced-by-markduplicates
    def __init__(self, *args, **kwargs):
        super(DuplicationHist, self).__init__(*args, **kwargs)
        self.kw = {
            'figure': {
                'plot_width': 400, 'plot_height': 400,
                'title': "Return of investment",
                'title_text_font_size': "10pt",
            },
            'renderer': {
                'groups': ["Sample"],
                'x': "BIN", 'y': "VALUE",
                'legend': 'Sample',
                'color': 'blue',
                'line_width': 2,
            },
            'xaxis': {'axis_label': 'Coverage multiple',
                      'axis_label_text_font_size': '10pt'},
            'yaxis': {'axis_label': 'Multiple of additional coverage',
                      'axis_label_text_font_size': '10pt'}}


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


def make_picard_summary_plots(inputfiles, ncol=4):
    d = {}
    TOOLS = "pan,box_zoom,wheel_zoom,box_select,lasso_select,resize,reset,save,hover"
    for (metrics_file, hist_file) in zip(inputfiles[0::2], inputfiles[1::2]):
        df_met = _read_metrics(metrics_file)
        df_hist = _read_metrics(hist_file)
        p1 = df_met.plot_metrics(tools=TOOLS)
        key = os.path.splitext(metrics_file)[0]
        if df_met.label not in d:
            d[df_met.label] = {}
            d[df_met.label][key] = {}
            d[df_met.label][key]['uri'] = [data_uri(metrics_file)]
            d[df_met.label][key]['file'] = [metrics_file]
        if df_hist is not None:
            p2 = df_hist.plot_hist(tools=TOOLS)
            d[df_met.label][key]['uri'].append(data_uri(hist_file))
            d[df_met.label][key]['file'].append(hist_file)
        else:
            p2 = []
        plist = p1 + p2
        gp = gridplot([plist[i:i+ncol] for i in range(0, len(plist), ncol)])
        d[df_met.label][key]['fig'] = gp
    return d
