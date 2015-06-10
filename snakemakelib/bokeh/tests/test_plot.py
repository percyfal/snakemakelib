# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904, C0301, C0103
import unittest
import pandas as pd
from nose.tools import raises
from bokeh.plotting import output_file, show
from bokeh.models import ColumnDataSource
from snakemakelib.bokeh.plot import make_dotplot, _data_fields


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
        make_dotplot(y='value')

    def test_dotplot_single_y_df_no_groups(self):
        p = make_dotplot(y='value', df=self.df, x='samples', groups=['_ALL'])
        (df, x, y) = _data_fields(p.children[0][0])
        #self.assertEqual(x, "_ALL")
        #self.assertEqual(y, "
        print (x, y)

        print( df)
        output_file("tabort.html")
        show(p)

    def test_dotplot_single_y_df_groups(self):
        p = make_dotplot(y='value', df=self.df, groups=['samples'])
        output_file("tabort.html")
        show(p)

    def test_dotplot_single_y_source_no_groups(self):
        p = make_dotplot(y='value', source=self.source)
        output_file("tabort.html")
        show(p)

    def test_dotplot_single_y_source_groups(self):
        p = make_dotplot(y='value', source=self.source, groups=['samples'])
        output_file("tabort.html")
        show(p)

    def test_dotplot_multiple_y_source_no_groups(self):
        p = make_dotplot(y=['value', 'count'], source=self.source)
        output_file("tabort.html")
        show(p)

    def test_dotplot_multiple_y_source_groups(self):
        p = make_dotplot(y=['value', 'count'], source=self.source,
                         groups=['samples'])
        output_file("tabort.html")
        show(p)

    def test_dotplot_multiple_y_df_groups_sidelegend(self):
        p = make_dotplot(y=['value', 'count'], source=self.source,
                         groups=['samples'], sidelegend=True)
        output_file("tabort.html")
        show(p)

    def test_dotplot_multiple_y_df_groups_sidelegend_relative(self):
        p = make_dotplot(y=['value', 'count'], source=self.source,
                         groups=['samples'], sidelegend=True,
                         relative_to="value")
        output_file("tabort.html")
        show(p)
