# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import pandas as pd
from bokeh.plotting import output_file
from bokeh.models import ColumnDataSource
from snakemakelib.bokeh.plot import make_dotplot


class TestPlot(unittest.TestCase):
    """Test bokeh plotting wrappers"""
    def setUp(self):
        self.df = pd.DataFrame([['foo', 0.23, 1, 'control'],
                                ['bar', 1.21, 3, 'case']],
                               columns=['samples', 'value',
                                        'count', 'treatment'])
        self.source = ColumnDataSource(self.df)

    def test_dotplot_single_y(self):
        #p = make_dotplot(y='value', df=self.df)
        #output_file("tabort.html")
        #show(p)
        print(dir(self.source))
        df = self.source.to_df()
        print(self.df)
        print(df)
        grouped = df.groupby('treatment')
        for name, group in grouped:
            print (name, group)
