# Copyright (C) 2015 by Per Unneberg
import math
import inspect
import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource
from bokeh.models.renderers import GlyphRenderer
from bokeh.palettes import brewer
from bokeh.plotting import figure, Figure, gridplot
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)


# Move elsewhere?
class QCArgs(object):
    """QC arguments"""

    def __init__(self, x, y, plot_type="line", **kwargs):
        self._x = x
        self._y = y
        self._plot_type = plot_type
        self._kwargs = kwargs

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def plot_type(self):
        return self._plot_type

    @property
    def kwargs(self):
        return self._kwargs


def _data_fields(p):
    """data_fields - get the df,x,y fields used by GlyphRenderer"""
    df = None
    for r in p.renderers:
        print(type(r))
        if isinstance(r, GlyphRenderer):
            spec = r.vm_serialize()
            df = spec['data_source'].to_df()
            spec = r.glyph.vm_serialize()
            x = spec['x']['field']
            y = spec['y']['field']
            break
    if df is None:
        return (None, None, None)
    return (df[[x, y]], x, y)


def _abline(p, slope=0, intercept=0, **kwargs):
    """abline - add an abline to current plot"""
    (df, x, y) = _data_fields(p)
    x0 = 0
    y0 = intercept
    x1 = max(df[x])
    y1 = (x1-x0) * slope + y0
    kwargs['color'] = kwargs.get('color', 'red')
    p.line(x=[x0, x1], y=[y0, y1], **kwargs)


def _sidelegend(y, color, **kwargs):
    df_leg = pd.DataFrame({'i': 1,
                           'y': list(range(min(len(y), len(color)))),
                           'text': y})
    source = ColumnDataSource(df_leg)
    lgd = figure(x_range=[0.5, 2.5], y_range=[-.5, len(color)-.5],
                 plot_width=kwargs.get('plot_width', 200),
                 plot_height=kwargs.get('plot_height', 400),
                 title="Legend",
                 title_text_font_size=kwargs.get('title_text_font_size',
                                                 "12"),
                 x_axis_type=None, y_axis_type=None)
    lgd.square(x=1, y='y', size=20, color=color, source=source)
    lgd.text(x=1.3, y='y', text='text', source=source, text_font_size="12")
    lgd.grid.grid_line_color = None
    return lgd


def _prepare_data(df=None, source=None):
    """Prepare data for plotting.

    Motivation: when grouping on several levels we need an extra index
    to plot on, otherwise plots will be off one index. Also, we check
    for existence of data frame and source, raising an error if none
    are provided.

    Args:
      df (py:class:`~pd.DataFrame`): a pandas data frame
      source (py:class:`~bokeh.models.ColumnDataSource`): a bokeh column
        data source
      groups (list): list of group names

    Returns:
      df, source: updated df and source
    """
    if source is None and df is None:
        raise TypeError(__name__ + """: both source and df None;
        you need to pass at least a data source or a data frame""")

    def _add_i():
        if '_i' in df.keys():
            raise KeyError
        if '_ALL' in df.keys():
            raise KeyError
        df['_i'] = list(range(1, len(df.index) + 1))
        df['_ALL'] = "ALL"

    if source is not None:
        df = source.to_df()
        # This is actually problematic; source should be immutable,
        # otherwise linked brushing won't work; however we need the
        # "_ALL" key for grouping. This points to the importance of
        # preparing the data *outside* the plots and not having any
        # defaults for x or groups
        _add_i()
    elif source is None:
        _add_i()
        source = ColumnDataSource(df)
    else:
        raise
    return df, source


# Model specialized plotting functions on bokeh.plotting.figure for simplicity
def scatterplot(x, y,
                xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                yaxis={'axis_label': "", 'major_label_orientation': 1},
                grid={'grid_line_color': None, 'grid_line_alpha': 1.0},
                tooltips=[], qc={}, **kwargs):

    """Make a scatter plot"""
    # x_axis_type, y_axis_type missed by inspect
    figmembers = [m[0] for m in inspect.getmembers(Figure)]
    + ['x_axis_type', 'y_axis_type']
    circle_kwargs_keys = list(set(list(kwargs.keys())).difference(set(figmembers)))
    circle_kwargs = {k: kwargs.pop(k) for k in circle_kwargs_keys}
    # Get a reference to a figure
    fig = figure(**kwargs)
    # Add qc cutoff if required
    if qc is not None:
        getattr(fig, qc.plot_type)(x=qc.x, y=qc.y, **qc.kwargs)
    # Add scatter points - check if list of strings, if present in
    # source we do several
    if set(y) <= set(circle_kwargs['source'].column_names):
        color = ["black"]
        if 'color' in circle_kwargs.keys():
            color = circle_kwargs.pop('color')
        for (yy, c) in zip(y, color):
            fig.circle(x=x, y=yy, color=c, legend=yy, **circle_kwargs)
    else:
        fig.circle(x=x, y=y, **circle_kwargs)
    # Apparently fig.xaxis and friends are lists so need to loop here.
    for k in xaxis.keys():
        [setattr(attr, k, xaxis[k]) for attr in fig.xaxis]
    for k in yaxis.keys():
        [setattr(attr, k, yaxis[k]) for attr in fig.yaxis]
    for k in grid.keys():
        [setattr(attr, k, grid[k]) for attr in fig.grid]
    # NB: currently assume it is a dictionary
    for tt in tooltips:
        h = fig.select(dict(type=tt['type']))
        h.tooltips = tt['tips']
    return (fig)


def make_dotplot(y, x="_i", df=None, source=None, groups=[], both=False,
                 xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                 yaxis={'axis_label': "", 'major_label_orientation': 1},
                 grid={'grid_line_color': "gray", 'grid_line_alpha': 0.3},
                 circle={'line_color': 'black', 'size': 8}, relative_to=None,
                 tooltips=[], sidelegend=False, **kwargs):
    """Make a (categorical) dotplot.

    Args:
      y (list): column name for y variable
      x (str): column name for x variable
      df (py:class:`~pandas.DataFrame`): a pandas data frame
      source (py:class:`~bokeh.models.ColumnDataSource`): column data source
      groups (list): list of column names to group by
      xaxis (dict): args passed to figure xaxis object
      yaxis (dict): args passed to figure yaxis object
      grid (dict): args passed to figure grid object
      circle (dict): args passed to figure circle object
      tooltips (list): tooltips to include
      relative_to (str): if set, calculate y values relative to parameter
      both (bool): if relative_to, return both plots
      sidelegend (bool): place legend beside figure
      kwargs (dict): keyword arguments to pass to figure method

    Returns:
      fig (py:class:`~bokeh.models.plots.GridPlot`): gridplot object
    """
    def _make_plot():
        color = brewer["PiYG"][min(max(3, len(y)), 10)]
        if 'color' in circle.keys():
            color = circle.pop('color')
        fig = figure(**kwargs)
        for (yy, c) in zip(y, color):
            if not sidelegend:
                circle.update({'legend': yy})
            fig.circle(x=x, y=yy, color=c,
                       source=source, **circle)
        for k in xaxis.keys():
            [setattr(attr, k, xaxis[k]) for attr in fig.xaxis]
        for k in yaxis.keys():
            [setattr(attr, k, yaxis[k]) for attr in fig.yaxis]
        for k in grid.keys():
            [setattr(attr, k, grid[k]) for attr in fig.grid]
        # NB: currently assume it is a dictionary
        for tt in tooltips:
            h = fig.select(dict(type=tt['type']))
            h.tooltips = tt['tips']
        return fig, color

    plist = []
    df, source = _prepare_data(df, source)
    grouped = df.groupby(groups) if groups else df.groupby('_ALL')
    if len(groups) > 1:
        kwargs['x_range'] = ["_".join(k)
                             for k in sorted(list(grouped.groups.keys()))]
    else:
        kwargs['x_range'] = list(grouped.groups.keys())
    # if not groups:
    #     x = "_ALL"
    if isinstance(y, str):
        y = [y]
    p, color = _make_plot()
    plist = [p]
    if relative_to is not None:
        df[y] = 100.0 * df[y].div(df[relative_to], axis="index")
        source = ColumnDataSource(df)
        kwargs['y_range'] = [0, max(110, max(df[y].max()))]
        yaxis['axis_label'] = "Proportion of {} (%)".format(relative_to)
        circle.update({'x_range': p.x_range})
        prel, color = _make_plot()
        prel.x_range = p.x_range
        plist.append(prel)
    # Make legend
    if sidelegend:
        lgd = _sidelegend(y, color)
        plist.append(lgd)
    return gridplot([plist])


def make_gridplot(y, df=None, source=None, x="_i",
                  text={'text': None, 'text_font_size': "6"}, groups=[],
                  ncol=1, total_plot_width=1200,
                  xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                  yaxis={'axis_label': "", 'major_label_orientation': 1},
                  grid={'grid_line_color': "gray", 'grid_line_alpha': 0.2},
                  circle={'line_color': 'black'}, line={},
                  share_x_range=False, share_y_range=False,
                  tooltips=[], qc=None, abline={}, **kwargs):
    """Make a gridplot.

    Args:
      y (str): column name for y variable
      df (py:class:`~pandas.DataFrame`): data frame
      source(py:class:`~bokeh.models.ColumnDataSource`): bokeh source table
      x (str): column name for x variable
      text (dict): args passed to figure text object,
                   with column name to use as text marker
      groups (list): list of column names to group by
      ncol(int): number of columns to group by
      total_plot_width (int): total plot width for a grid line
      xaxis (dict): args passed to figure xaxis object
      yaxis (dict): args passed to figure yaxis object
      grid (dict): args passed to figure grid object
      share_x_range (bool): share x range
      share_y_range (bool): share y range
      tooltips (list): tooltips to include
      qc (bool): add QC indicators
      abline (dict): add abline to plot, where parameters to abline
                      must be intercept and slope, and possibly color
      kwargs (dict): keyword arguments to pass to figure method

    Returns:
      gp (py:class:`bokeh.models.GridPlot`): GridPlot object

    """
    (df, source) = _prepare_data(df, source)
    first = True
    plist = []
    grouped = df.groupby(groups) if groups else df.groupby("_ALL")
    # if not groups:
    #     x = "_ALL"
    kwargs['plot_width'] = int(total_plot_width / ncol)
    kwargs['plot_height'] = int(total_plot_width / ncol)
    title = [kwargs.pop('title')] if 'title' in kwargs else []
    for name, data in grouped:
        # NB: sources cannot and should not be shared here!
        source = ColumnDataSource(data)
        p = figure(title=", ".join(title + [name]),
                   **kwargs)
        if text['text'] is not None:
            p.text(x=x, y=y, source=source, **text)
        elif line:
            p.line(x=x, y=y, source=source, **line)
        else:
            p.circle(x=x, y=y, source=source, **circle)
        if abline:
            _abline(p, **abline)
        for k in xaxis.keys():
            [setattr(attr, k, xaxis[k]) for attr in p.xaxis]
        for k in yaxis.keys():
            [setattr(attr, k, yaxis[k]) for attr in p.yaxis]
        for k in grid.keys():
            [setattr(attr, k, grid[k]) for attr in p.grid]
        plist.append(p)
        # Update kwargs for subsequent plots
        if share_x_range and first:
            kwargs['x_range'] = plist[0].x_range
        if share_y_range and first:
            kwargs['y_range'] = plist[0].y_range
        first = False
    gp = gridplot([plist[i:i+ncol] for i in range(0, len(plist), ncol)])
    return gp


def make_scatterplot(y, df=None, source=None, x="_i", groups=[],
                     xaxis={'axis_label': "",
                            'major_label_orientation': np.pi/3},
                     yaxis={'axis_label': "",
                            'major_label_orientation': 1},
                     grid={'grid_line_color': "gray", 'grid_line_alpha': 0.2},
                     circle={'line_color': 'black'}, tooltips=[],
                     **kwargs):
    """Make a scatter plot

    Args:
      y (str): column name for y variable
      df (py:class:`~pandas.DataFrame`): data frame
      source(py:class:`~bokeh.models.ColumnDataSource`): bokeh source table
      x (str): column name for x variable
      groups (list): list of column names to group by
      xaxis (dict): args passed to figure xaxis object
      yaxis (dict): args passed to figure yaxis object
      grid (dict): args passed to figure grid object
      circle (dict): args passed to figure circle object
      tooltips (list): tooltips to include
      kwargs (dict): keyword arguments to pass to figure method

    Returns:
      fig (py:class:`bokeh.models.Figure`): Figure object
    """
    df, source = _prepare_data(df, source)
    grouped = df.groupby(groups) if groups else df.groupby("_ALL")
    colors = list({k: v for (k, v) in
                   zip(grouped.groups.keys(), brewer["PiYG"]
                       [min(max(3, len(grouped.groups.keys())), 10)]
                       * grouped.size()[0])}.values()) * grouped.size()[0]
    # Fix ranges
    x_range = list(df[x]) if df[x].dtype is np.dtype('object') else []
    if 'x_range' in kwargs:
        x_range = kwargs.pop('x_range')

    # Get figure
    fig = figure(x_range=x_range, **kwargs)
    # Add scatter points - check if list of strings, if present in
    # source we do several
    if set(y) <= set(source.column_names):
        color = brewer["PiYG"][10]
        if 'color' in circle.keys():
            color = circle.pop('color')
        for (yy, c) in zip(y, color):
            fig.circle(x=x, y=yy, color=c, legend=yy,
                       source=source, **circle)
    else:
        fig.circle(x=x, y=y, source=source, color=colors,
                   legend=groups, **circle)
    for k in xaxis.keys():
        [setattr(attr, k, xaxis[k]) for attr in fig.xaxis]
    for k in yaxis.keys():
        [setattr(attr, k, yaxis[k]) for attr in fig.yaxis]
    for k in grid.keys():
        [setattr(attr, k, grid[k]) for attr in fig.grid]
    # NB: currently assume it is a dictionary
    for tt in tooltips:
        h = fig.select(dict(type=tt['type']))
        h.tooltips = tt['tips']
    return fig


def make_lineplot(df=None, source=None, x=None, y=None, groups=[],
                  xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                  yaxis={'axis_label': "", 'major_label_orientation': 1},
                  grid={'grid_line_color': None, 'grid_line_alpha': 1.0},
                  line={}, tooltips=[], qc={}, **kwargs):
    """Make a scatter plot

    Args:
      df (py:class:`~pandas.DataFrame`): data frame
      source(py:class:`~bokeh.models.ColumnDataSource`): bokeh source table
      x (str): column name for x variable
      y (str): column name for y variable
      groups (list): list of column names to group by
      xaxis (dict): args passed to figure xaxis object
      yaxis (dict): args passed to figure yaxis object
      grid (dict): args passed to figure grid object
      line (dict): args passed to figure line object
      tooltips (list): tooltips to include
      qc (dict): qc parameters
      kwargs (dict): keyword arguments to pass to figure method

    Returns:
      fig (py:class:`bokeh.models.Figure`): Figure object
    """
    df, source = _prepare_data(df, source)
    fig = figure(**kwargs)
    grouped = df.groupby(groups)
    colors = {k: v for (k, v) in
              zip(grouped.groups.keys(), brewer["PiYG"]
                  [min(max(3, len(grouped.groups.keys())), 10)] *
                  math.ceil(len(grouped.groups.keys()) / 10))}
    for i in grouped.groups.keys():
        labels = grouped.get_group(i).columns
        xname = labels[0]
        # Here we currently assume second column is y; this is not
        # always the case for insertion metrics
        yname = labels[1]
        x = getattr(grouped.get_group(i), xname)
        y = getattr(grouped.get_group(i), yname)
        fig.line(x, y, legend=i, color=colors[i], **line)
    for k in xaxis.keys():
        [setattr(attr, k, xaxis[k]) for attr in fig.xaxis]
    for k in yaxis.keys():
        [setattr(attr, k, yaxis[k]) for attr in fig.yaxis]
    for k in grid.keys():
        [setattr(attr, k, grid[k]) for attr in fig.grid]
    return fig
