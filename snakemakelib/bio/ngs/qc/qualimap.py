# Copyright (C) 2015 by Per Unneberg
import pandas as pd
import numpy as np
from snakemake.report import data_uri
from snakemakelib.results import Results
from snakemakelib.bokeh.plot import make_gridplot, make_dotplot
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
        df_tmp['chrlen_percent'] = 100 * df_tmp['chrlen'] /\
            sum(df_tmp['chrlen'])
        df_tmp['mapped_bases_percent'] = 100 * df_tmp['mapped_bases'] /\
            sum(df_tmp['mapped_bases'])
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


def make_qualimap_plots(qmglobals=None, coverage_per_contig=None,
                        **kwargs):
    """Make qualimap summary plots"""
    retval = {'fig': {'coverage_per_contig': None, 'globals': None},
              'file': {'coverage_per_contig': coverage_per_contig,
                       'globals': qmglobals},
              'uri': {'coverage_per_contig': data_uri(coverage_per_contig),
                      'globals': data_uri(qmglobals)}}
    # Globals
    if qmglobals is not None:
        df_all = pd.read_csv(qmglobals)
        df_all['number of unique reads'] = df_all['number of mapped reads']\
            - df_all['number of duplicated reads']
        plot_config = {'y': ['number of reads',
                             'number of mapped reads',
                             'number of duplicated reads',
                             'number of unique reads'],
                       'df': df_all, 'groups': ['Sample'],
                       'y_range': [0, max(df_all['number of reads'])],
                       'xaxis': {'axis_label': "sample",
                                 'major_label_orientation': np.pi/3,
                                 'axis_label_text_font_size': '10pt'},
                       'yaxis': {'axis_label': "count",
                                 'major_label_orientation': 1,
                                 'axis_label_text_font_size': '10pt'},
                       'circle': {'size': 20, 'alpha': 0.3},
                       'title': "Qualimap summary",
                       'relative_to': "number of reads",
                       'plot_width': 400, 'plot_height': 400,
                       'both': True}
        retval['fig']['globals'] = make_dotplot(**plot_config)

    # Coverage per contig
    if coverage_per_contig is not None:
        df_all = pd.read_csv(coverage_per_contig, index_col=0)
        gp = make_gridplot(x="chrlen_percent", y="mapped_bases_percent",
                           df=df_all, groups=['Sample'], ncol=4,
                           share_x_range=True, share_y_range=True,
                           xaxis={'axis_label': "Chromosome length of total (%)",
                                  'axis_label_text_font_size': "8"},
                           yaxis={'axis_label': "Mapped bases of total (%)",
                                  'axis_label_text_font_size': "8"},
                           kwtext={'text': 'chr', 'text_font_size': "2"},
                           abline={'slope': 1, 'intercept': -1, 'pad': 1.5},
                           title_text_font_size="10", **kwargs)
        retval['fig']['coverage_per_contig'] = gp

    return retval
