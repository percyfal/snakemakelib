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
                else:
                    df.append(df_tmp, ignore_index=True)
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
    source = ColumnDataSource(df_all)
    df = df_all
    p = figure(title="plot",
               x_range=[-1, max(df['chrlen_percent']) * 1.1],
               y_range=[-1, max(df['mapped_bases_percent']) * 1.1])
    # p.text('chrlen_percent',
    #        'mapped_bases_percent',
    #        text='chr', source=source)
    p.circle('chrlen_percent',
             'mapped_bases_percent',
             source=source)
    # p.line(x=[0, max('chrlen_percent')],
    #        y=[0, max('mapped_bases_percent')],
    #        source=source)
    callback = Callback(args=dict(source=source), code="""
    var dsource = source.get('data');
    var ind = cb_obj.get('value')
    var arrayUnique = function(a) {
        return a.reduce(function(p, c) {
            if (p.indexOf(c) < 0) p.push(c);
                return p;
            }, []);
    };
    var unique_samples = arrayUnique(dsource['Sample'])
    var sample = unique_samples[ind - 1]
    samples = dsource['Sample']
    x = dsource['x']
    y = dsource['y']
    for (i = 0; i < samples.length; i++) {
         if (sample == samples[i]) {
    x[i] = 2* x[i]
         }
    }
    source.trigger('change');
    """)
    # dropdown = Dropdown(menu=[(s, "change") for s in samples],
    #                     callback=callback, label="Samples")
    slider = Slider(start=0, end=len(samples), value=0, step=1,
                    title="Sample", callback=callback)
    return {'fig': vform(slider, p),
            'uri': data_uri(coverage_per_contig),
            'file': coverage_per_contig}
