# Copyright (C) 2015 by Per Unneberg
import os
import pandas as pd
import jinja2
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, BoxSelectTool
from bokeh.models.widgets import VBox, HBox, TableColumn, DataTable
from bokeh.plotting import figure, output_file, show, gridplot
from bokeh.palettes import brewer
from snakemake.report import data_uri
from snakemakelib.report.utils import recast, trim_header
from snakemakelib.bokeh.plot import scatterplot, QCArgs, scatterplot2

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

def make_star_alignment_plots(inputfile, do_qc=False, min_reads=200000, min_map=40, max_unmap=20):
    """Make star alignment plots"""
    df = pd.read_csv(inputfile, index_col=0)
    samples = list(df.index)
    # Currently hover tool and categorical variables don't play
    # nicely together in bokeh: see
    # https://github.com/bokeh/bokeh/issues/624

    # Workaround as long as categorical variables don't work with HoverTool
    df['i'] = list(range(0, len(df.index)))
    df['samples'] = samples
    df['mismatch_sum'] = df['Mismatch_rate_per_base__PCT'] + df['Deletion_rate_per_base'] + df['Insertion_rate_per_base']
    df['PCT_of_reads_unmapped'] = df['PCT_of_reads_unmapped:_other'] + df['PCT_of_reads_unmapped:_too_many_mismatches'] + df['PCT_of_reads_unmapped:_too_short']

    colors = brewer["PiYG"][3]
    colormap = {'False' : colors[0], 'True' : colors[1]}
    
    columns = [
        TableColumn(field="samples", title="Sample"),
        TableColumn(field="Number_of_input_reads", title="Number of input reads"),
        TableColumn(field="Uniquely_mapped_reads_PCT", title="Uniquely mapped reads (%)"),
        TableColumn(field="Mismatch_rate_per_base__PCT", title="Mismatch rate per base (%)"),
        TableColumn(field="Insertion_rate_per_base", title="Insertion rate per base (%)"),
        TableColumn(field="Deletion_rate_per_base", title="Deletion rate per base (%)"),
        TableColumn(field="PCT_of_reads_unmapped", title="Unmapped reads (%)"),
    ]
        
    source = ColumnDataSource(df)
    # Generate the table
    table = DataTable(source=source, columns=columns, editable=False, width = 1000)

    # Default tools, plot_config and tooltips
    TOOLS="pan,wheel_zoom,box_zoom,box_select,lasso_select,reset,save,hover"
    plot_config=dict(plot_width=400, plot_height=400, tools=TOOLS, title_text_font_size='12pt',
                     x_axis_type = 'linear', x_range = [0, len(samples)],
                     xaxis = {'axis_label' : 'sample', 'axis_label_text_font_size' : '10pt', 'major_label_orientation' : np.pi/3},
                     yaxis = {'axis_label' : 'reads', 'axis_label_text_font_size' : '10pt', 'major_label_orientation' : np.pi/3}
                     )

    # Number of input reads
    c1 = list(map(lambda x: colormap[str(x)], df['Number_of_input_reads'] < min_reads)) if do_qc else "blue"
    qc = QCArgs(x=[0,len(samples)], y=[min_reads, min_reads], line_dash=[2,4]) if do_qc else None
    p1 = scatterplot(x='i', y='Number_of_input_reads', source=source, color=c1, qc=qc,
                     title="Number of input reads",
                     tooltips = [{'type':HoverTool, 'tips' : [('Sample', '@samples'),('Reads', '@Number_of_input_reads'),]}],
                     y_range=[0, max(df['Number_of_input_reads'])],
                     y_axis_type = "log",
                     **plot_config)

    # Uniquely mapped reads
    plot_config.update({'y_axis_type' : 'linear', 'axis_label' : 'percent (%)'})
    c2 = list(map(lambda x: colormap[str(x)], df['Uniquely_mapped_reads_PCT'] < min_map))  if do_qc else "blue"
    qc = QCArgs(x=[0,len(samples)], y=[min_map, min_map], line_dash=[2,4]) if do_qc else None
    p2 = scatterplot(x='i', y='Uniquely_mapped_reads_PCT', source=source, color=c2, qc=qc,
                     title="Uniquely mapping reads",
                     y_range = [0, 100],
                     tooltips = [{'type':HoverTool, 'tips' : [('Sample', '@samples'),('Pct_mapped', '@Uniquely_mapped_reads_PCT'),]}],
                     **plot_config)

    # Mapping reads in general
    c3 = list(map(lambda x: colormap[str(x)], df['PCT_of_reads_unmapped'] > max_unmap))  if do_qc else "blue"
    qc = QCArgs(x=[0,len(samples)], y=[max_unmap, max_unmap], line_dash=[2,4]) if do_qc else None
    p3 = scatterplot(x='i', y='PCT_of_reads_unmapped',
                     source=source, color=c3, qc=qc, title="Unmapped reads",
                     y_range = [0, 100],
                     tooltips = [{'type':HoverTool, 'tips' : [('Sample', '@samples'),('Pct_unmapped', '@PCT_of_reads_unmapped'),]}], **plot_config)
    
    # Mismatch/indel rate
    plot_config['tools'] = TOOLS.replace("lasso_select,", "")
    plot_config['yaxis'].update({'axis_label' : 'Rate per base'})
    p4 = scatterplot(x='i', y = ['Mismatch_rate_per_base__PCT', 'Insertion_rate_per_base', 'Deletion_rate_per_base'],
                color = ["blue", "red", "green"], source = source,
                title = "Mismatch and indel rates",
                tooltips =  [{'type':HoverTool, 'tips' : [('Sample', '@samples'),
                                                            ('Mismatch rate per base', '@Mismatch_rate_per_base__PCT'),
                                                            ('Insertion rate per base', '@Insertion_rate_per_base'),
                                                            ('Deletion rate per base', '@Deletion_rate_per_base'),
                                                          ]},
                                                          ],
                **plot_config)
    select_tool = p4.select(dict(type=BoxSelectTool))
    select_tool.dimensions=['width']

    # Plot sum
    plot_config['yaxis'].update({'axis_label' : 'Mismatch/indel sum'})
    c5 = list(map(lambda x: colormap[str(x)], df['mismatch_sum'] > 1.0))  if do_qc else "blue"
    qc = QCArgs(x=[0,len(samples)], y=[1.0, 1.0], line_dash=[2,4]) if do_qc else None
    p5 = scatterplot(x='i', y='mismatch_sum',
                     source=source, color=c5, qc=qc, title="Mismatch / indel sum",
                     tooltips = [{'type':HoverTool, 'tips' : [('Sample', '@samples'),('Mismatch/indel rate per base', '@mismatch_sum'),]}], **plot_config)
    select_tool = p5.select(dict(type=BoxSelectTool))
    select_tool.dimensions=['width']

    # Plot histogram of ratio
    # plot_config['tools'] = "pan,box_zoom,reset,save"
    # p6 = figure(title="Histogram of mismatch and indel rates", **plot_config)
    # hist, edges = np.histogram(df['mismatch_sum'], density=False, bins=50)
    # p6.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
    #    fill_color="#036564", line_color="#033649")
    # p6.xaxis.axis_label = "Mismatch/indel sum"
    # p6.yaxis.axis_label = "Count"

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
    
    return {'fig' : VBox(children=[gridplot([[p1, p2, p3]]), HBox(children=[gridplot([[p4, p5]])])]),
            'table' : table, 'qctable' : df_qc,
            'uri' : data_uri(inputfile),
            'file' : inputfile}
