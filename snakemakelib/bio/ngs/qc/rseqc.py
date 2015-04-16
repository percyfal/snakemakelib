# Copyright (C) 2015 by Per Unneberg
import os
import re
import io
import pandas as pd
import jinja2
import numpy as np
from bokeh.models import HoverTool, ColumnDataSource, BoxSelectTool
from bokeh.models.widgets import VBox, HBox, TableColumn, DataTable
from bokeh.plotting import figure, output_file, show, gridplot
from bokeh.palettes import brewer
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

def collect_rseqc_results(inputfiles, samples):
    """Collect rseqc results"""
    first = True
    df_rd = df_gc = None
    for (f, s) in zip(inputfiles, samples):
        smllogger.debug("parsing input file {f} for sample {s}".format(f=f, s=s))
        try:
            infile = os.path.join(os.path.dirname(f), "read_distribution.txt")
            with open(infile, "r") as fh:
                data = fh.readlines()
            d1 = io.StringIO("".join([re.sub(r"(\w) (\w)", r"\1_\2", x) for x in data[3:6]]))
            tmp  = pd.read_table(d1, engine="python", sep="[ ]+", header=None)
            tmp.columns = ["Group", "Tag_count"]
            tmp["Total_bases"] = "NA"
            tmp["Tags/Kb"] = "NA"
            d2 = io.StringIO("".join(data[7:17]))
            df_rd_tmp = pd.read_table(d2, engine="python", sep="[ ]+")
            df_rd_tmp = df_rd_tmp.append(tmp)
            df_rd_tmp["sample"] = s
        except:
            smllogger.warn("pandas reading table {f} failed".format(f=infile))
        try:
            infile = os.path.join(os.path.dirname(f), "geneBody_coverage.geneBodyCoverage.txt")
            df_gc_tmp = pd.read_table(infile, engine="python", sep="\t")
            df_gc_tmp.index = [s]
        except:
            smllogger.warn("pandas reading table {f} failed".format(f=infile))
        try:
            if first:
                df_rd = df_rd_tmp
                df_gc = df_gc_tmp
                first = False
            else:
                df_rd = df_rd.append(df_rd_tmp)
                df_gc = df_gc.append(df_gc_tmp)
        except:
            smllogger.warn("failed to append data")
    return {'rd': df_rd, 'gc':df_gc}

def make_rseqc_summary_plots(df_rd, df_gc, do_qc=True, min_exonmap=60.0, max_three_prime_map=10.0):
    """Make rseqc summary plots"""
    samples = list(df_gc.index)
    # Use tags for formula 
    df = df_rd.pivot_table(columns=["Group"], values=["Tag_count"], index="sample")
    df['Tag_count', "ExonMap"] = 100.0 * (df['Tag_count', "CDS_Exons"] + df['Tag_count', "3'UTR_Exons"] + df['Tag_count', "5'UTR_Exons"]) / df['Tag_count', "Total_Assigned_Tags"]

    df.columns = df.columns.droplevel()
    df['i'] = list(range(0, len(df.index)))
    df['samples'] = samples
    df_gc["three_prime_map"] = 100.0 * df_gc.loc[:, "91":"100"].sum(axis=1) / df_gc.loc[:, "1":"100"].sum(axis=1)
    df = pd.concat([df, df_gc], axis=1)

    colors = brewer["PiYG"][3]
    colormap = {'False' : colors[0], 'True' : colors[2]}
    columns = [
        TableColumn(field="samples", title="Sample"),
        TableColumn(field="ExonMap", title="Tags mapping to exons (%)"),
        TableColumn(field="3' Map", title="Tags mapping to 3' end (%)"),
    ]
    source = ColumnDataSource(df)

    # Default tools, plot_config and tooltips
    TOOLS="pan,box_zoom,box_select,lasso_select,reset,save,hover"
    plot_config=dict(plot_width=300, plot_height=300, tools=TOOLS, title_text_font_size='12pt')
    
    # Exonmap plot
    c1 = list(map(lambda x: colormap[str(x)], df['ExonMap'] < min_exonmap)) if do_qc else colors[0]
    p1 = figure(x_range=[0, len(samples)], y_range=[0, 105], x_axis_type=None, y_axis_type="linear", title="Tags mapping to exons", **plot_config)
    if do_qc:
        p1.line(x=[0,len(samples)], y=[min_exonmap, min_exonmap], line_dash=[2,4])
    p1.circle(x='i', y='ExonMap', color=c1, source=source)
    p1.xaxis.major_label_orientation = np.pi/3
    p1.grid.grid_line_color = None
    hover = p1.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ('ExonMap', '@ExonMap'),
    ]

    # Fraction reads mapping to the 10% right-most end
    c2 = list(map(lambda x: colormap[str(x)], df['three_prime_map'] > max_three_prime_map)) if do_qc else colors[0]
    p2 = figure(x_range=[0, len(samples)], y_range=[0, 105], x_axis_type=None, y_axis_type="linear", title="Reads mapping to 3' end", **plot_config)
    p2.circle(x='i', y='three_prime_map', color=c2, source=source)
    if do_qc:
        p2.line(x=[0,len(samples)], y=[max_three_prime_map, max_three_prime_map], line_dash=[2,4])
    p2.xaxis.major_label_orientation = np.pi/3
    p2.grid.grid_line_color = None
    hover = p2.select(dict(type=HoverTool))
    hover.tooltips = [
        ('Sample', '@samples'),
        ("3' map", '@three_prime_map'),
    ]

    return {'plots' : gridplot([[p1, p2]])}
