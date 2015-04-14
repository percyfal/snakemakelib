# Copyright (C) 2015 by Per Unneberg
import os
import pandas as pd
import jinja2
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, BoxSelectTool
from bokeh.models.widgets import VBox, TableColumn, DataTable
from bokeh.plotting import figure, output_file, show, gridplot
from snakemakelib.report.utils import recast, trim_header

def collect_star_alignment_results(input, samples):
    """Collect star alignment results"""
    df = None
    for (f, s) in zip(input, samples):
        df_tmp = pd.read_table(f, sep="|", names=["name", "value"], engine="python", skiprows=[7,22,27])
        d = {trim_header(x, underscore=True, percent=True):recast(y) for (x,y) in zip(df_tmp["name"], df_tmp["value"])}
        if df is None:
            df = pd.DataFrame(d, index=[s])
        else:
            df = df.append(pd.DataFrame(d, index=[s]))
    return df

def make_star_alignment_plots(df, samples, min_reads=200000, min_map=40):
    """Make star alignment plots"""
    # Currently hover tool and categorical variables don't play
    # nicely together in bokeh: see
    # https://github.com/bokeh/bokeh/issues/624

    # Workaround as long as categorical variables don't work with HoverTool
    df['i'] = list(range(0, len(df.index)))
    df['samples'] = samples
    
    # Subset data frame for cutoffs
    df['reads_filtered'] = df['Number_of_input_reads'] < min_reads
    df['map_filtered'] = df['Uniquely_mapped_reads_PCT'] < min_map
    df['filtered'] = (df['Uniquely_mapped_reads_PCT'] < min_map) | (df['Number_of_input_reads'] < min_reads)

    colormap = {'False':'blue', 'True':'red'}
    
    columns = [
        TableColumn(field="samples", title="Sample"),
        TableColumn(field="Number_of_input_reads", title="Number of input reads"),
        TableColumn(field="Uniquely_mapped_reads_PCT", title="Uniquely mapped reads, PCT"),
    ]
        
    source = ColumnDataSource(df)
    # Generate the table
    table = DataTable(source=source, columns=columns, editable=False, width = 1000)

    # Default tools, plot_config and tooltips
    TOOLS="pan,box_zoom,box_select,lasso_select,reset,save,hover"
    plot_config=dict(plot_width=300, plot_height=300, tools=TOOLS, title_text_font_size='10pt')

    # Number of input reads
    c1 = list(map(lambda x: colormap[str(x)], df['reads_filtered']))
    p1 = figure(x_range=[0, len(samples)], x_axis_type=None, y_axis_type="log", title="Number of input reads", **plot_config)
    p1.line(x=[0,len(samples)], y=[min_reads, min_reads], line_dash=[2,4])
    p1.circle(x='i', y='Number_of_input_reads', color=c1, source=source)
    p1.xaxis.major_label_orientation = np.pi/3
    p1.grid.grid_line_color = None
    hover = p1.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ('Num_input_reads', '@Number_of_input_reads'),
    ]

    # Uniquely mapped reads
    c2 = list(map(lambda x: colormap[str(x)], df['map_filtered']))
    p2 = figure(x_range=[0, len(samples)], y_range=[-5,105], x_axis_type=None, title="Uniquely mapping reads", **plot_config)
    p2.line(y=[min_map,min_map], x=[0, len(samples)], line_dash=[2,4])
    p2.circle(x='i', y='Uniquely_mapped_reads_PCT', color=c2, source=source)
    p2.xaxis.major_label_orientation = np.pi/3
    p2.grid.grid_line_color = None
    hover = p2.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ('Pct_mapped_reads', '@Uniquely_mapped_reads_PCT'),
    ]


    # Uniquely mapped reads vs number of input reads
    c3 = list(map(lambda x: colormap[str(x)], df['filtered']))
    p3 = figure(title="Number of input reads vs uniquely mapping reads", y_axis_type="log", **plot_config)
    p3.circle(y='Number_of_input_reads', x='Uniquely_mapped_reads_PCT', color=c3, source=source)
    p3.xaxis.major_label_orientation = np.pi/3
    p3.grid.grid_line_color = None
    hover = p3.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ('Num_input_reads', '@Number_of_input_reads'),
        ('Pct_mapped_reads', '@Uniquely_mapped_reads_PCT'),
    ]

    return {'plots' : gridplot([[p1, p2, p3]]), 'table' : table}
