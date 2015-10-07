# Copyright (C) 2015 by Per Unneberg
import re
import os
import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from Bio import SeqIO
from bokeh.plotting import figure, gridplot
from bokeh.io import vform
from bokeh.models import ColumnDataSource, DataTable, TableColumn, HoverTool, BoxSelectTool, CustomJS
from bokeh.models.widgets import Select, Toggle
from bokehutils.axes import xaxis, yaxis
from bokehutils.geom import dotplot, points
from bokehutils.tools import tooltips
from bokehutils.color import colorbrewer
from snakemakelib.log import LoggerManager

logger = LoggerManager().getLogger(__name__)

TOOLS = "pan,wheel_zoom,box_zoom,box_select,lasso_select,resize,reset,save,hover"

def scrnaseq_alignment_qc_plots(rseqc_read_distribution=None, rseqc_gene_coverage=None,
                      star_results=None):
    """Make alignment QC plots for scrnaseq workflow

    Args:
      rseqc_read_distribution (str): RSeQC read distribution results csv file
      rseqc_gene_coverage (str): RSeQC gene coverage results csv file
      star_results (str): star alignment results csv file

    Returns: 
      dict: dictionary with keys 'fig' pointing to a (:py:class:`~bokeh.models.GridPlot`) Bokeh GridPlot object and key 'table' pointing to a (:py:class:`~bokeh.widgets.DataTable`) DataTable

    """
    df_star = pd.read_csv(star_results, index_col="Sample")
    df_rseqc_rd = pd.read_csv(rseqc_read_distribution, index_col="Sample").reset_index().pivot_table(columns=["Group"], values=["Tag_count"], index=["Sample"])
    df_rseqc_rd.columns = ["_".join(x) if isinstance(x, tuple) else x for x in df_rseqc_rd.columns]
    df_rseqc_gc = pd.read_csv(rseqc_gene_coverage, index_col="Sample")
    df_all = df_star.join(df_rseqc_rd)
    df_all = df_all.join(df_rseqc_gc['three_prime_map'])
    source = ColumnDataSource(df_all)
    columns = [
        TableColumn(field="Sample", title="Sample"),
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
    kwfig = {'plot_width': 400, 'plot_height': 400, 
             'title_text_font_size': "12pt"}
    kwxaxis = {'axis_label': 'Sample',
               'major_label_orientation': np.pi/3}
    kwyaxis = {'axis_label_text_font_size': '10pt',
               'major_label_orientation': np.pi/3}

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

    # Fraction reads mapping to 10% right-most end
    p7 = figure(title="Reads mapping to 3' end",
                x_range=p1.x_range,
                tools=TOOLS,
                **kwfig)
    dotplot(p7, "Sample", "three_prime_map", source=source)
    xaxis(p7, **kwxaxis)
    yaxis(p7, axis_label="Percent", **kwyaxis)
    tooltips(p7, HoverTool, [('Sample', '@Sample'),
                             ("3' map", '@three_prime_map'), ])


    return {'fig': gridplot([[p1, p2, p3], [p4, p5, p6], [p7, None, None]]),
            'table': table}

# FIXME: move to stats module
def pca(expr, **kwargs):
    """scrnaseq pca - run pca

    Args:
      expr (DataFrame): pandas data frame with expression values, rows
      correspond to samples/experimental units, columns to genes
      kwargs: keyword  arguments

    Returns:
      pcaobj (PCA): PCA model fitted to expr
    """
    pcaobj = PCA(n_components=kwargs.get("n_components", 10))
    pcaobj.fit(expr)
    return pcaobj

# FIXME: move to stats module
def pca_results(pcaobj, expr, metadata=None, **kwargs):
    """Generate output file on which to base pca plots
    
    Args:
      pcaobj (PCA): pca model
      expr (DataFrame): pandas data frame with expression values
      metadata (DataFrame): sample metadata

    Returns:
      pcares (DataFrame): pca results concatenated,
      loadings (DataFrame): pca loadings of top 10 units for each component
    """
    pcares = pd.DataFrame(pcaobj.fit(expr).transform(expr))
    if not expr.index.name is None:
        pcares.index = expr.index
    if not metadata is None:
        md = pd.read_csv(metadata, index_col=0)
        pcares = pcares.join(md)
    # Fix for bokehutils; all columns names must be strings...
    pcares.columns = [str(x) for x in list(pcares.columns)]
    return pcares

def scrnaseq_pca_plots(pca_results_file=None, metadata=None, pcaobjfile=None):
    """Make PCA QC plots for scrnaseq workflow

    Args:
      pca_results_file (str): pca results file
      metadata (str): metadata file name
      pcaobjfile (str): file name containing pickled pca object

    Returns: 
      dict: dictionary with keys 'fig' pointing to a (:py:class:`~bokeh.models.GridPlot`) Bokeh GridPlot object and key 'table' pointing to a (:py:class:`~bokeh.widgets.DataTable`) DataTable

    """
    if not metadata is None:
        md = pd.read_csv(metadata, index_col=0)
    if not pcaobjfile is None:
        with open(pcaobjfile, 'rb') as fh:
            pcaobj = pickle.load(fh)
    df_pca = pd.read_csv(pca_results_file, index_col="sample")
    df_pca['color'] = ['red'] * df_pca.shape[0]
    df_pca['x'] = df_pca['0']
    df_pca['y'] = df_pca['1']
    df_pca['size'] = [10] * df_pca.shape[0]

    # TODO: Move to stats.pca
    # Define the loadings
    # loadings = pd.DataFrame(pcaobj.components_).T
    # sort the columns
    #    loadings_index = 
    # Annotate the genes
    

    source = ColumnDataSource(df_pca)

    # Add radio button group
    cmap = colorbrewer(datalen = df_pca.shape[0])
    callback  = CustomJS(args=dict(source=source), code="""
var data = source.get('data');
var active = cb_obj.get('active');
// If radio button group need index
// var label = cb_obj.get('label')[active]
var label = cb_obj.get('label');
var Paired = {
    3	: [ "#A6CEE3","#1F78B4","#B2DF8A"],
    4	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C"],
    5	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99"],
    6	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C"],
    7	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F"],
    8	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00"] ,
    9	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6"],
    10	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A"],
    11	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99"],
    12	: [ "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99", "#B15928"]};
var colormap = {};
// Comment out if(!active) if radiobuttongroup
if (!active) {
    if (typeof(data[label][0]) != "string") {
	var minsize = Math.log(Math.min.apply(null, data[label]))/Math.LN2;
	var maxsize = Math.log(Math.max.apply(null, data[label]))/Math.LN2;
	for (i = 0; i < data['sample'].length; i++) {
	    data['size'][i] = 10 + (Math.log(data[label][i])/Math.LN2 - minsize) / (maxsize - minsize) * 20;
	}
    } else {
	var j = 0;
	// Simple categorical grouping
	for (i = 0; i < data['sample'].length; i++) {
            if (data[label][i] in colormap) {
            } else {
		colormap[data[label][i]] = j;
		j++;
            }
	}
	var nfac = Object.keys(colormap).length;
	if (nfac > 12) {
            nfac = 12;
	} 
	if (nfac < 3) {
	    nfac = 3;
	}
	var colors = Paired[nfac];
	for (i = 0; i < data[label].length; i++) {
            data['color'][i] = colors[colormap[data[label][i]]]
	}
    }
} else {
    if (typeof(data[label][i] != "string")) {
	for (i = 0; i < data['sample'].length; i++) {
	    data['size'][i] = 10;
	}
    }
}
source.trigger('change');
    """)
    if not md is None:
        # Waiting for callbacks to be implemented upstream in bokeh
        # rbg = RadioButtonGroup(labels=list(md.columns),
        #                        callback=callback)
        toggle_buttons = [Toggle(label=x, callback=callback) for x in list(md.columns) + ["TPM", "FPKM"]]
    else:
        toggle_buttons = []
        # rbg = RadioButtonGroup()
    # PC components
    xcallback = CustomJS(args=dict(source=source), code="""
    var myRe = /^([0-9]+)/;
        var data = source.get('data');
        var value = myRe.exec(String(cb_obj.get('value')))[1]
        x = data['x']
        for (i = 0; i < x.length; i++) {
              x[i] = data[value][i]
        }
        source.trigger('change');
    """)
    ycallback = CustomJS(args=dict(source=source), code="""
    var myRe = /^([0-9]+)/;
        var data = source.get('data');
        var value = myRe.exec(String(cb_obj.get('value')))[1]
        y = data['y']
        for (i = 0; i < y.length; i++) {
             y[i] = data[value][i]
        }
        source.trigger('change');
    """)

    pca_components = sorted([int(x) + 1 for x in source.column_names if re.match("\d+", x)])
    menulist = ["{} ({:.2f}%)".format(x, 100.0 * p) for x, p in zip(pca_components, pcaobj.explained_variance_ratio_)]
    component_x = Select(title = "PCA component x", options = menulist, value=menulist[0],
                         callback=xcallback)
    component_y = Select(title = "PCA component y", options = menulist, value=menulist[1],
                         callback=ycallback)

    buttons = toggle_buttons + [component_x, component_y]
    # Make the pca plot
    kwfig = {'plot_width': 400, 'plot_height': 400,
             'title_text_font_size': "12pt"}


    p1 = figure(title="Principal component analysis",
                tools=TOOLS, **kwfig)

    points(p1, 'x', 'y', source=source, color='color', size='size',
           alpha=.7)
    kwxaxis = {
               'axis_label_text_font_size': '10pt',
               'major_label_orientation': np.pi/3}
    kwyaxis = {
        'axis_label_text_font_size': '10pt',
        'major_label_orientation': np.pi/3}
    xaxis(p1, **kwxaxis)
    yaxis(p1, **kwyaxis)
    tooltiplist = [("sample", "@sample")] if "sample" in source.column_names else []
    if not md is None:
        tooltiplist = tooltiplist + [(str(x), "@{}".format(x)) for x
                                     in md.columns] + \
        [("Detected genes (TPM)", "@TPM"), ("Detected genes (FPKM)", "@FPKM")]
    tooltips(p1, HoverTool, tooltiplist)

    # Detected genes, FPKM and TPM
    p2 = figure(title="Number of detected genes",
                x_range=list(df_pca.index), tools=TOOLS,
                **kwfig)
    kwxaxis.update({'axis_label': "Sample"})
    kwyaxis.update({'axis_label': "Detected genes"})
    dotplot(p2, "sample", "FPKM", source=source)
    xaxis(p2, **kwxaxis)
    yaxis(p2, **kwyaxis)
    tooltips(p2, HoverTool, [('sample', '@sample'),
                             ('# genes (FPKM)', '@FPKM')])
    return {'pca' : vform(*(buttons + [gridplot([[p1, p2]])]))}


def scrnaseq_extra_ref_plot(expr=None, extra_ref=None, unit_id="gene_id"):
    """Plot spikeins.

    FIXME: this should be a generic plot for multi-plotting gene
    expression values over samples

    Args: 
      expr (str): expression matrix file name (long format)
      extra_ref (list): list of extra reference file names
      unit_id (str): name of column to use for plotting

    """
    
    if extra_ref is None:
        return
    # Add spikein plots
    extra_ids = []
    seqio_format = {'fa' : 'fasta',
                    'fasta' : 'fasta'}
    for ref in extra_ref:
        fh = open(ref)
        ext = os.path.splitext(ref)[1].replace(".", "")
        extra_ids += [record.id.strip() for record in SeqIO.parse(fh, format=seqio_format[ext])]
    expr_long = pd.read_csv(expr)
    criterion = expr_long["gene_id"].map(lambda x: x in extra_ids)
    import math
    expr_long.loc[criterion, "TPM_log"] = [math.log2(x + 1) for x in expr_long.loc[criterion, "TPM"]]
    expr_long.loc[criterion, "size"] = 5
    df = expr_long.loc[criterion, :]
    fig = figure(title="spikein plot", x_range=sorted(list(set(df["sample"]))),
                 plot_width=800, plot_height=400, tools=TOOLS)
    cmap = colorbrewer(size=4, datalen = len(list(set(df["sample"]))))
    for (groupid, g), c in zip(df.groupby('gene_id'), cmap):
        source = ColumnDataSource(g)
        fig.circle(x='sample', y='TPM_log', legend=groupid, color=c,
                   source=source, size=g['size'])
    xaxis(fig, axis_label='Sample', major_label_orientation= np.pi/3)
    yaxis(fig, axis_label='TPM_log', major_label_orientation= np.pi/3)
    tooltips(fig, HoverTool,[('sample', '@sample'),
                             (unit_id, '@' + unit_id)])
    # Widget for selecting ref
    # selectref_callback = CustomJS(args=dict(source=source), code="""
    #     var data = source.get('data');
    #     var value = cb_obj.get('value')
    #     y = data['y']
    #     for (i = 0; i < y.length; i++) {
    #          y[i] = data[value][i]
    #     }
    #     }

    #     source.trigger('change');
    # """)
    # selectref = Select(title = "Reference", options = extra_ids,
    #                    value=extra_ids[0], callback=selectref_callback)

    #return {'spikein' : vform(*([selectref] + [fig]))}
    return {'spikein' : fig}

def scrnaseq_expr_heatmap(expr_file=None, unit_id="gene_id", metadata=None):
    from sklearn.cluster import AgglomerativeClustering
    if not expr_file:
        return None
    expr_long = pd.read_csv(expr_file)
    expr = expr_long.pivot_table(columns=unit_id, values="TPM",
                                 index="sample")
    import math
    expr_log = expr.applymap(lambda x: math.log2(x + 1))
    n_samples = expr_log.shape[0]
    samplecorr = expr_log.T.corr()
    model = AgglomerativeClustering(linkage='average',
                                    affinity='precomputed',
                                    n_clusters=n_samples)

    fit = model.fit(samplecorr)
    order = model.fit_predict(samplecorr)
    samplecorr_reordered = samplecorr.iloc[order, order]
    print(samplecorr_reordered.head())
    sc_stack = samplecorr_reordered.stack().reset_index()
    print (sc_stack.head())
    sc_stack.index = samplecorr_reordered["sample"]
    sc_stack.columns = ['x', 'y', 'cor']
    print (md.head())
    print (samplecorr_reordered.head())
    if not metadata is None:
        md = pd.read_csv(metadata, index_col=0)
        samplecorr_reordered = samplecorr_reordered.join(md)
    print(samplecorr_reordered.head())
    # from bokeh.charts import Scatter
    # fig = Scatter(expr
    # print(dir(model))
    

    
