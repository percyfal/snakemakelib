# Copyright (C) 2015 by Per Unneberg
import pandas as pd
import numpy as np
from bokeh.plotting import figure, output_file, show, gridplot
from bokeh.models import ColumnDataSource, DataTable, TableColumn, HoverTool, BoxSelectTool
from bokehutils.axes import xaxis, yaxis
from bokehutils.geom import dotplot
from bokehutils.tools import tooltips

def scrnaseq_qc_plots(rseqc_read_distribution=None, rseqc_gene_coverage=None,
                      star_results=None):
    """Make QC plots for scrnaseq workflow

    Args:
      rseqc_read_distribution (str): RSeQC read distribution results csv file
      rseqc_gene_coverage (str): RSeQC gene coverage results csv file
      star_results (str): star alignment results csv file

    Returns:
      gp (:py:class:`~bokeh.models.GridPlot`): Bokeh GridPlot object
    """
    df_star = pd.read_csv(star_results, index_col="Sample")
    df_rseqc_rd = pd.read_csv(rseqc_read_distribution, index_col="Sample").reset_index().pivot_table(columns=["Group"], values=["Tag_count"], index=["Sample"])
    df_rseqc_rd.columns = ["_".join(x) if isinstance(x, tuple) else x for x in df_rseqc_rd.columns]
    df_rseqc_gc = pd.read_csv(rseqc_gene_coverage, index_col="Sample")
    df_all = df_star.join(df_rseqc_rd)
    source = ColumnDataSource(df_all)
    columns = [
        TableColumn(field="samples", title="Sample"),
        TableColumn(field="Number_of_input_reads",
                    title="Number of input reads"),
        TableColumn(field="Uniquely_mapped_reads_PCT",
                    title="Uniquely mapped reads (%)"),
        TableColumn(field="Mismatch_rate_per_base__PCT",
                    title="Mismatch rate per base (%)"),
        TableColumn(field="Insertion_rate_per_base",
                    title="Insertion rate per base (%)"),
        TableColumn(field="Deletion_rate_per_base",
                    title="Deletion rate per base (%)"),
        TableColumn(field="PCT_of_reads_unmapped",
                    title="Unmapped reads (%)"),
    ]
    table = DataTable(source=source, columns=columns,
                      editable=False, width=1000)
    TOOLS = "pan,wheel_zoom,box_zoom,box_select,lasso_select,resize,reset,save,hover"
    kwfig = {'plot_width': 400, 'plot_height': 400, 
             'title_text_font_size': "12pt"}
    kwxaxis = {'axis_label': 'Sample',
               'major_label_orientation': np.pi/3}
    kwyaxis = {'axis_label_text_font_size': '10pt',
               'major_label_orientation': np.pi/3}
    output_file("tabort.html")

    # Input reads
    p1 = figure(title="Number of input reads",
                x_range=list(df_all.index), tools=TOOLS,
                y_axis_type="log", **kwfig)
    dotplot(p1, "Sample", "Number_of_input_reads", source=source)
    xaxis(p1, **kwxaxis)
    yaxis(p1, axis_label="Reads", **kwyaxis)
    tooltips(p1, HoverTool, [('Sample', '@Sample'),
                             ('Reads', '@Number_of_input_reads')])

    # Uniquely mapping
    p2 = figure(title="Uniquely mapping reads",
                x_range=p1.x_range,
                y_range=[0, 100],
                tools=TOOLS,
                **kwfig)
    dotplot(p2, "Sample", "Uniquely_mapped_reads_PCT", source=source)
    xaxis(p2, **kwxaxis)
    yaxis(p2, axis_label="Percent", **kwyaxis)
    tooltips(p2, HoverTool, [('Sample', '@Sample'),
                             ('Pct_mapped', '@Uniquely_mapped_reads_PCT')])

    # Unmapped
    p3 = figure(title="Unmapped reads",
                x_range=p1.x_range,
                y_range=[0, 100],
                tools=TOOLS,
                **kwfig)
    dotplot(p3, "Sample", "PCT_of_reads_unmapped", source=source)
    xaxis(p3, **kwxaxis)
    yaxis(p3, axis_label="Percent", **kwyaxis)
    tooltips(p3, HoverTool, [('Sample', '@Sample'),
                             ('Pct_unmapped', '@PCT_of_reads_unmapped')])

    # Mismatch/indel rate
    p4 = figure(title="Mismatch and indel rates",
                x_range=p1.x_range,
                tools=TOOLS,
                **kwfig)
    dotplot(p4, "Sample", "Mismatch_rate_per_base__PCT", source=source, legend="Mismatch")
    dotplot(p4, "Sample", "Insertion_rate_per_base", source=source, legend="Insertion", color="red")
    dotplot(p4, "Sample", "Deletion_rate_per_base", source=source, legend="Deletion", color="green")
    xaxis(p4, **kwxaxis)
    yaxis(p4, axis_label="Percent", **kwyaxis)
    tooltips(p4, HoverTool,  [('Sample', '@samples'),
                                              ('Mismatch rate per base',
                                               '@Mismatch_rate_per_base__PCT'),
                                              ('Insertion rate per base',
                                               '@Insertion_rate_per_base'),
                                              ('Deletion rate per base',
                                               '@Deletion_rate_per_base'), ])
    select_tool = p4.select(dict(type=BoxSelectTool))
    select_tool.dimensions = ['width']

    

             
    # Unmapped
    p5 = figure(title="Mismatch/indel sum",
                x_range=p1.x_range,
                tools=TOOLS,
                **kwfig)
    dotplot(p5, "Sample", "mismatch_sum", source=source)
    xaxis(p5, **kwxaxis)
    yaxis(p5, axis_label="Percent", **kwyaxis)
    tooltips(p5, HoverTool, [('Sample', '@Sample'),
                             ('Mismatch/indel rate per base',
                              '@mismatch_sum'), ])
    select_tool = p5.select(dict(type=BoxSelectTool))
    select_tool.dimensions = ['width']

    # Fraction reads mapping to 10% right-most end
    p6 = figure(title="Tags mapping to exons",
                x_range=p1.x_range,
                tools=TOOLS,
                **kwfig)
    dotplot(p6, "Sample", "Tag_count_ExonMap", source=source)
    xaxis(p6, **kwxaxis)
    yaxis(p6, axis_label="Percent", **kwyaxis)
    tooltips(p6, HoverTool, [('Sample', '@Sample'),
                             ('ExonMap', '@Tag_count_ExonMap'), ])


    show(gridplot([[p1, p2, p3], [p4, p5, p6]]))
             
