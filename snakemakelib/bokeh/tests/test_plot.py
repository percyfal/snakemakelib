# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103, C0111, R0201, E1101
import unittest
import pandas as pd
from nose.tools import raises
from bokeh.models import ColumnDataSource
from bokeh.plotting import output_file, show
from snakemakelib.bokeh.plot import make_dotplot, _data_fields, make_gridplot
from snakemakelib.bokeh.plot import make_lineplot, make_scatterplot


class TestPlot(unittest.TestCase):
    """Test bokeh plotting wrappers"""
    def setUp(self):
        self.df = pd.DataFrame([['foo', 0.23, 1, 'control'],
                                ['bar', 1.21, 3, 'case']],
                               columns=['samples', 'value',
                                        'count', 'treatment'])
        self.source = ColumnDataSource(self.df)

    @raises(TypeError)
    def test_dotplot_no_source(self):
        make_dotplot(x='samples', y='value')

    def test_dotplot_single_y_df_no_groups(self):
        p = make_dotplot(y='value', df=self.df, x='samples')
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(y, ['value'])
        self.assertTupleEqual(df.shape, (2, 2))

    def test_dotplot_single_y_df_groups(self):
        p = make_dotplot(x='samples', y='value', df=self.df, groups=['samples'])
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(y, ['value'])
        self.assertTupleEqual(df.shape, (2, 2))

    def test_dotplot_single_y_source_no_groups(self):
        p = make_dotplot(x='samples', y='value', source=self.source)
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(y, ['value'])
        self.assertTupleEqual(df.shape, (2, 2))

    def test_dotplot_single_y_source_groups(self):
        p = make_dotplot(x='samples', y='value', source=self.source, groups=['samples'])
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(y, ['value'])
        self.assertTupleEqual(df.shape, (2, 2))

    def test_dotplot_multiple_y_source_no_groups(self):
        p = make_dotplot(x='samples', y=['value', 'count'], source=self.source)
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(sorted(y), sorted(['value', 'count']))
        self.assertTupleEqual(df.shape, (2, 3))

    def test_dotplot_multiple_y_source_groups(self):
        p = make_dotplot(x='samples', y=['value', 'count'],
                         source=self.source,
                         groups=['samples'])
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(sorted(y), sorted(['value', 'count']))
        self.assertTupleEqual(df.shape, (2, 3))

    def test_dotplot_multiple_y_df_groups_sidelegend(self):
        p = make_dotplot(x='samples', y=['value', 'count'],
                         source=self.source, groups=['samples'],
                         sidelegend=True)
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(sorted(y), sorted(['value', 'count']))
        self.assertTupleEqual(df.shape, (2, 3))
        self.assertEqual(len(p.children[0]), 2)

    def test_dotplot_multiple_y_df_groups_sidelegend_relative(self):
        p = make_dotplot(x='samples', y=['value', 'count'],
                         source=self.source, groups=['samples'],
                         sidelegend=True, relative_to="value")
        (df, x, y) = _data_fields(p.children[0][0])
        self.assertListEqual(x, ['samples'])
        self.assertListEqual(sorted(y), sorted(['value', 'count']))
        self.assertTupleEqual(df.shape, (2, 3))
        self.assertEqual(len(p.children[0]), 3)

    def test_gridplot(self):
        p = make_gridplot(x='samples', y='value', df=self.df,
                          groups=['treatment'], plot_width=400,
                          plot_height=400, ncol=3)
        output_file("tabort.html")
        show(p)

    def test_scatterplot(self):
        p = make_scatterplot(x='samples', y=['count', 'value'], df=self.df,
                             groups=['treatment'], plot_width=400,
                             plot_height=400)
        output_file("tabort.html")
        show(p)

    def test_lineplot(self):
        p = make_lineplot(x='count', y='value', df=self.df,
                          groups=['samples'], plot_width=400,
                          plot_height=400)
        output_file("tabort.html")
        show(p)
        
