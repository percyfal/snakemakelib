# Copyright (C) 2015 by Per Unneberg
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokehutils.geom import points, abline
from bokehutils.mgeom import mdotplot
from bokehutils.facet import facet_grid
from bokehutils.axes import xaxis, yaxis, main
from snakemake.report import data_uri
from snakemakelib.results import Results
from snakemakelib.log import LoggerManager


smllogger = LoggerManager().getLogger(__name__)

COVERAGE_PER_CONTIG_COLUMNS = ["chr", "chrlen", "mapped_bases",
                               "mean_coverage", "sd"]
GLOBALS_COLUMNS = ["name", "value"]


class Qualimap(Results):
    _keys = ['globals', 'coverage_per_contig']

    def __init__(self, *args, **kw):
        super(Qualimap, self).__init__(*args, **kw)

    def _collect_globals(self, data, first, sample):
        df_tmp = self.parse_data(data,
                                 rs=("Globals", "Insert"),
                                 skip=1, split=True,
                                 columns=GLOBALS_COLUMNS,
                                 dtype=float, sep=" = ")
        df_tmp['value'] = [float(x.split(" ")[0].replace(",", ""))
                           for x in df_tmp['value']]
        df_tmp['Sample'] = sample
        try:
            if first:
                self['globals'] = df_tmp
            else:
                self['globals'] = self['globals'].append(df_tmp, ignore_index=True)
        except:
            smllogger.warn("failed to append data to globals dataframe")

    def _collect_coverage_per_contig(self, data, first, sample):
        df_tmp = self.parse_data(data,
                                 rs=("Coverage per contig", None),
                                 skip=2, split=True,
                                 columns=COVERAGE_PER_CONTIG_COLUMNS,
                                 dtype=float)
        df_tmp["Sample"] = sample
        try:
            df_tmp['chrlen_percent'] = 100 * df_tmp['chrlen'] /\
                sum(df_tmp['chrlen'])
            df_tmp['mapped_bases_percent'] = 100 * df_tmp['mapped_bases'] /\
                sum(df_tmp['mapped_bases'])
        except:
            smllogger.warn("coverage_per_contig: failed to normalize data")
        try:
            if first:
                self['coverage_per_contig'] = df_tmp
            else:
                self['coverage_per_contig'] = self['coverage_per_contig'].append(
                    df_tmp, ignore_index=True)
        except:
            smllogger.warn("failed to append data to coverage_per_contig dataframe")

    def _collect_results(self):
        smllogger.info("Collecting results")
        first = True
        for (f, s) in zip(self._inputfiles, self._samples):
            smllogger.debug("Reading input file {f} for sample {s}".format(f=f, s=s))
            data = self.load_lines(f)
            self._collect_globals(data, first, s)
            self._collect_coverage_per_contig(data, first, s)
            first = False
        if self['globals'] is not None:
            self['globals'] = self['globals'].pivot(
                index='Sample', columns='name', values='value')
            self['globals']['number of unique reads'] = self['globals']['number of mapped reads']\
                                                        - self['globals']['number of duplicated reads']



def make_qualimap_plots(qmglobals=None, coverage_per_contig=None):
    """Make qualimap summary plots"""
    retval = {'fig': {'coverage_per_contig': None, 'globals': None},
              'file': {'coverage_per_contig': coverage_per_contig,
                       'globals': qmglobals},
              'uri': {'coverage_per_contig': data_uri(coverage_per_contig),
                      'globals': data_uri(qmglobals)}}
    # Globals
    if qmglobals is not None:
        df_all = pd.read_csv(qmglobals)
        fig = figure(y_range=[0, max(df_all['number of reads'])],
                     title="Mapping summary",
                     title_text_font_size="12pt",
                     plot_width=400, plot_height=400,
                     x_range=list(set(df_all["Sample"])))
        mdotplot(fig, x='Sample', size=10,
                 df=df_all, alpha=0.5,
                 y= ['number of reads',
                     'number of mapped reads',
                     'number of duplicated reads',
                     'number of unique reads'])
        xaxis(fig, axis_label="sample",
              major_label_orientation=np.pi/3,
              axis_label_text_font_size='10pt')
        yaxis(fig, axis_label="count",
              major_label_orientation=1,
              axis_label_text_font_size='10pt')
        retval['fig']['globals'] = fig
        
    # Coverage per contig
    if coverage_per_contig is not None:
        df_all = pd.read_csv(coverage_per_contig, index_col=0)
        fig = figure(width=300, height=300)
        points(fig, x="chrlen_percent", y="mapped_bases_percent",
               df=df_all, glyph="text", text="chr", text_font_size="8pt")
        main(fig, title_text_font_size="8pt")
        xaxis(fig, axis_label="Chromosome length of total (%)",
              axis_label_text_font_size="8pt")
        yaxis(fig, axis_label="Mapped bases of total (%)",
              axis_label_text_font_size="8pt")

        gp = facet_grid(fig, x="chrlen_percent", y="mapped_bases_percent",
                        df=df_all, groups=["Sample"], width=300, height=300,
                        share_x_range=True, share_y_range=True,
                        title_text_font_size="12pt")
        for fig in [item for sublist in gp.children for item in sublist]:
            abline(fig, x="chrlen_percent", y="mapped_bases_percent", df=df_all, slope=1)
        retval['fig']['coverage_per_contig'] = gp
    return retval        
