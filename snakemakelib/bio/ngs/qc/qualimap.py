# Copyright (C) 2015 by Per Unneberg
import pandas as pd
from bokeh.io import vform
from bokeh.plotting import figure
from bokeh.models import Callback, Slider, ColumnDataSource, Dropdown
from snakemake.report import data_uri
from snakemakelib.results import Results
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

COVERAGE_PER_CONTIG_COLUMNS = ["chr", "chrlen", "mapped_bases",
                               "mean_coverage", "sd"]


class Qualimap(Results):
    _keys = ['coverage_per_contig']

    def __init__(self, *args, **kw):
        super(Qualimap, self).__init__(*args, **kw)

    def _collect_results(self):
        smllogger.info("Collecting results")
        first = True
        for (f, s) in zip(self._inputfiles, self._samples):
            print ("File: ", f)
            data = self.load_lines(f)
            df_tmp = self.parse_data(data,
                                     rs=("Coverage per contig", None),
                                     skip=2, split=True,
                                     columns=COVERAGE_PER_CONTIG_COLUMNS,
                                     dtype=float)
            df_tmp["Sample"] = s
            try:
                if first:
                    df = df_tmp.copy(deep=True)
                    first = False
                    print ("first: ", df.shape)
                else:
                    print (df.shape)
                    print (df.iloc[1,])
                    print (df_tmp.iloc[1,])
                    df = df.append(df_tmp, ignore_index=True)
                    print (df.shape)
            except:
                smllogger.warn("failed to append data to coverage_per_contig dataframe")
        df['chrlen_percent'] = 100 * df['chrlen']/sum(df['chrlen'])
        df['mapped_bases_percent'] = 100 * df['mapped_bases']/sum(df['mapped_bases'])
        self['coverage_per_contig'] = df


def make_qualimap_plots(coverage_per_contig=None,
                        **kwargs):
    """Make qualimap summary plots"""
    df_all = pd.read_csv(coverage_per_contig, index_col=0)
    samples = list(set(df_all.Sample))
    print (samples)
    print (list(df_all.chrlen_percent))
    plist = []
    source = ColumnDataSource(df_all)
    df = df_all
    p = figure(title="plot",
               x_range=[-1, max(df['chrlen_percent']) * 1.1],
               y_range=[-1, max(df['mapped_bases_percent']) * 1.1])
    p.text('chrlen_percent',
           'mapped_bases_percent',
           text='chr', source=source)
    m = max(df['chrlen_percent'] + df['mapped_bases_percent'])
    p.line(x=[0, m], y=[0, m])
    callback = Callback(args=dict(source=source), code="""
    var data = source.get('data');
    var sample = cb_obj.get('value')
    x = data['mapped_bases_percent']
    y = data['chrlen']
    text = data['chr']
    samples = data['Sample']
    for (i = 0; i < samples.length; i++) {
         if (sample == samples[i]) {
             x[i] = data['mapped_bases_percent'][i]
             y[i] = data['chrlen_percent'][i]
             text[i] = data['chr'][i]
         }
    }
    source.trigger('change');
    """)
    dropdown= Dropdown(menu = [(s, s) for s in samples], callback=callback, label="Samples")
    return {'fig': vform(dropdown, p),
            'uri': data_uri(coverage_per_contig),
            'file': coverage_per_contig}
