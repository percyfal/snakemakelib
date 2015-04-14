# Copyright (C) 2015 by Per Unneberg
import os
import pandas as pd
import jinja2
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, BoxSelectTool
from bokeh.models.widgets import VBox, HBox, TableColumn, DataTable
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

def make_star_alignment_plots(df, samples, do_qc=False, min_reads=200000, min_map=40):
    """Make star alignment plots"""
    # Currently hover tool and categorical variables don't play
    # nicely together in bokeh: see
    # https://github.com/bokeh/bokeh/issues/624

    # Workaround as long as categorical variables don't work with HoverTool
    df['i'] = list(range(0, len(df.index)))
    df['samples'] = samples
    df['mismatch_sum'] = df['Mismatch_rate_per_base__PCT'] + df['Deletion_rate_per_base'] + df['Insertion_rate_per_base']
    colormap = {'False':'blue', 'True':'red'}
    
    columns = [
        TableColumn(field="samples", title="Sample"),
        TableColumn(field="Number_of_input_reads", title="Number of input reads"),
        TableColumn(field="Uniquely_mapped_reads_PCT", title="Uniquely mapped reads (%)"),
        TableColumn(field="Mismatch_rate_per_base__PCT", title="Mismatch rate per base (%)"),
        TableColumn(field="Insertion_rate_per_base", title="Insertion rate per base (%)"),
        TableColumn(field="Deletion_rate_per_base", title="Deletion rate per base (%)"),
    ]
        
    source = ColumnDataSource(df)
    # Generate the table
    table = DataTable(source=source, columns=columns, editable=False, width = 1000)

    # Default tools, plot_config and tooltips
    TOOLS="pan,box_zoom,box_select,lasso_select,reset,save,hover"
    plot_config=dict(plot_width=300, plot_height=300, tools=TOOLS, title_text_font_size='12pt')
    # Number of input reads
    c1 = list(map(lambda x: colormap[str(x)], df['Number_of_input_reads'] < min_reads)) if do_qc else "blue"
    p1 = figure(x_range=[0, len(samples)], x_axis_type=None, y_axis_type="log", title="Number of input reads", **plot_config)
    if do_qc:
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
    c2 = list(map(lambda x: colormap[str(x)], df['Uniquely_mapped_reads_PCT'] < min_map))  if do_qc else "blue"
    p2 = figure(x_range=[0, len(samples)], y_range=[-5,105], x_axis_type=None, title="Uniquely mapping reads", **plot_config)
    if do_qc:
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
    c3 = list(map(lambda x: colormap[str(x)], (df['Uniquely_mapped_reads_PCT'] < min_map) | (df['Number_of_input_reads'] < min_reads)))  if do_qc else "blue"
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

    # Mismatch/indel rate
    plot_config['tools'] = TOOLS.replace("lasso_select,", "")
    p4 = figure(title="Mismatch and indel rates", x_axis_type=None, **plot_config)
    p4.circle(x='i', y='Mismatch_rate_per_base__PCT', color="blue", source=source)
    p4.circle(x='i', y='Insertion_rate_per_base', color="red", source=source)
    p4.circle(x='i', y='Deletion_rate_per_base', color="green", source=source)
    p4.xaxis.major_label_orientation = np.pi/3
    p4.grid.grid_line_color = None
    hover = p4.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ('Mismatch rate per base', '@Mismatch_rate_per_base__PCT'),
        ('Insertion rate per base', '@Insertion_rate_per_base'),
        ('Deletion rate per base', '@Deletion_rate_per_base'),
    ]
    select_tool = p4.select(dict(type=BoxSelectTool))
    select_tool.dimensions=['width']
    p4.yaxis.axis_label = "Rate per base"
    p4.xaxis.axis_label = "Sample"

    # Plot sum
    c5 = list(map(lambda x: colormap[str(x)], df['mismatch_sum'] > 1.0))  if do_qc else "blue"
    p5 = figure(title="Mismatch/indel sum", **plot_config)
    p5.circle(x='i', y='mismatch_sum', color=c5, source=source)
    select_tool = p5.select(dict(type=BoxSelectTool))
    select_tool.dimensions=['width']
    hover = p5.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ('Mismatch/indel rate per base', '@mismatch_sum'),
    ]
    p5.yaxis.axis_label = "Mismatch/indel sum"

    # Plot histogram of ratio
    plot_config['tools'] = "pan,box_zoom,reset,save"
    p6 = figure(title="Histogram of mismatch and indel rates", **plot_config)
    hist, edges = np.histogram(df['mismatch_sum'], density=False, bins=50)
    p6.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
       fill_color="#036564", line_color="#033649")
    p6.xaxis.axis_label = "Mismatch/indel sum"
    p6.yaxis.axis_label = "Count"

    df_qc = None
    if do_qc:
        # QC summary table
        d = {'samples':samples,
            'read_filter' : df['Number_of_input_reads'] < min_reads,
            'map_filter' : df['Uniquely_mapped_reads_PCT'] < min_map,
            'mismatch_filter' : df['mismatch_sum'] > 1.0,
            }
        d['filter'] = d['read_filter'] | d['map_filter'] | d['mismatch_filter']
        df_qc = pd.DataFrame(data=d, index=df.samples)
    
    return {'plots' : VBox(children=[gridplot([[p1, p2, p3]]), HBox(children=[gridplot([[p4, p5]]), p6])]), 'table' : table, 'qctable' : df_qc}
