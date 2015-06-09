# Copyright (C) 2015 by Per Unneberg
import math
import inspect
import numpy as np
import pandas as pd
from bokeh.io import hplot
from bokeh.models import ColumnDataSource
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


def scatterplot2(y, df, x="i", groups=[],
                 xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                 yaxis={'axis_label': "", 'major_label_orientation': 1},
                 grid={'grid_line_color': None, 'grid_line_alpha': 1.0},
                 circle={}, tooltips=[], qc=None, **kwargs):
    """Make a scatter plot

    Args:
      y (str): column name for y variable
      df (py:class:`~pandas.DataFrame`): data frame
      x (str): column name for x variable
      groups (list): list of column names to group by
      xaxis (dict): args passed to figure xaxis object
      yaxis (dict): args passed to figure yaxis object
      grid (dict): args passed to figure grid object
      circle (dict): args passed to figure circle object
      tooltips (list): tooltips to include
      qc (bool): do qc
      kwargs (dict): keyword arguments to pass to figure method

    Returns:
      fig (py:class:`bokeh.models.Figure`): Figure object
    """
    # Update data frame
    df['i'] = list(range(1, len(df.index) + 1))
    # Catchall groupby group
    df['all'] = "ALL"
    # grouping
    g = df.groupby(groups) if groups else df.groupby('all')
    colors = list({k: v for (k, v) in
                   zip(g.groups.keys(), brewer["PiYG"]
                       [min(max(3, len(g.groups.keys())), 10)]
                       * g.size()[0])}.values()) * g.size()[0]

    source = ColumnDataSource(df)
    # Fix ranges
    x_range = list(df[x]) if df[x].dtype is np.dtype('object') else []
    x = "i" if df[x].dtype is np.dtype('object') else x

    # Get figure
    fig = figure(x_range=x_range, **kwargs)

    # Add qc cutoff if required
    if qc:
        getattr(fig, qc.plot_type)(x=qc.x, y=qc.y, **qc.kwargs)
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


def lineplot(df, x=None, y=None, groups=[],
             xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
             yaxis={'axis_label': "", 'major_label_orientation': 1},
             grid={'grid_line_color': None, 'grid_line_alpha': 1.0},
             tooltips=[], qc=None, **kwargs):
    fig = figure(**kwargs)
    g = df.groupby(groups)
    colors = {k: v for (k, v) in
              zip(g.groups.keys(), brewer["PiYG"]
                  [min(max(3, len(g.groups.keys())), 10)] *
                  math.ceil(len(g.groups.keys()) / 10))}
    for i in g.groups.keys():
        labels = g.get_group(i).columns
        xname = labels[0]
        # Here we currently assume second column is y; this is not
        # always the case for insertion metrics
        yname = labels[1]
        x = getattr(g.get_group(i), xname)
        y = getattr(g.get_group(i), yname)
        fig.line(x, y, legend=i, color=colors[i])
    for k in xaxis.keys():
        [setattr(attr, k, xaxis[k]) for attr in fig.xaxis]
    for k in yaxis.keys():
        [setattr(attr, k, yaxis[k]) for attr in fig.yaxis]
    return (fig)


def dotplot(y, df, groups=[],
            xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
            yaxis={'axis_label': "", 'major_label_orientation': 1},
            grid={'grid_line_color': None, 'grid_line_alpha': 1.0},
            tooltips=[], **kwargs):
    """Make a dotplot"""
    g = df.groupby(groups)
    df['i'] = list(range(1, len(df.index) + 1))
    source = ColumnDataSource(df.sort())
    figmembers = [x[0] for x in inspect.getmembers(Figure)]
    circle_kwargs_keys = list(set(list(kwargs.keys())).difference(
        set(figmembers)))
    circle_kwargs = {k: kwargs.pop(k) for k in circle_kwargs_keys}
    x_range = ["_".join(x) for x in sorted(list(g.groups.keys()))]
    fig = figure(x_range=x_range, **kwargs)
    # Applying jitter fails for some reason: j = -
    # math.ceil(len(y)/2); offset = 1 / (3*len(y))
    if set(y) <= set(source.column_names):
        color = ['black'] * len(y)
        color = brewer["PiYG"][10]  # [min(max(3, len(y)), 10)]
        if 'color' in circle_kwargs.keys():
            color = circle_kwargs.pop('color')
        for (yy, c) in zip(y, color):
            fig.circle(x="i", y=yy, color=c, legend=yy,
                       source=source, alpha=1, **circle_kwargs)
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


def make_dotplot(y, df, groups=[], both=False,
                 xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                 yaxis={'axis_label': "", 'major_label_orientation': 1},
                 grid={'grid_line_color': "gray", 'grid_line_alpha': 0.2},
                 circle={}, relative_to=None,
                 tooltips=[], sidelegend=False, **kwargs):
    """Make a (categorical) dotplot.

    Args:
      y (str): column name for y variable
      df (py:class:`~pandas.DataFrame`): data frame
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
      fig (py:class:`bokeh.models.Figure`): Figure object
    """
    def _make_plot():
        df['i'] = list(range(1, len(df.index) + 1))
        source = ColumnDataSource(df)
        color = brewer["PiYG"][min(max(3, len(y)), 10)]
        if 'color' in circle.keys():
            color = circle.pop('color')
        fig = figure(x_range=x_range, y_range=y_range, **kwargs)
        for (yy, c) in zip(y, color):
            if not sidelegend:
                circle.update({'legend': yy})
            fig.circle(x="i", y=yy, color=c,
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
    y_range = kwargs.pop('y_range')
    if groups:
        grouped = df.groupby(groups)
    else:
        grouped = {"all": df}
    if len(groups) > 1:
        x_range = ["_".join(x) for x in sorted(list(grouped.groups.keys()))]
    else:
        x_range = sorted(list(grouped.groups.keys()))
    p, color = _make_plot()
    plist = [p]
    if relative_to is not None:
        df_tmp = df[y].T
        iloc = next((i for i in range(len(y)) if relative_to in y[i]))
        df_tmp = (100.0 * df_tmp/df_tmp.iloc[iloc]).T
        df = df_tmp
        y_range = [0, 110]
        yaxis['axis_label'] = "Proportion of {} (%)".format(relative_to)
        circle.update({'x_range': p.x_range})
        prel, color = _make_plot()
        plist.append(prel)
    # Make legend
    if sidelegend:
        df_leg = pd.DataFrame({'i': 1,
                               'y': list(range(min(len(y), len(color)))),
                               'text': y})
        source = ColumnDataSource(df_leg)
        lgd = figure(x_range=[0.5, 2.5], y_range=[-.5, len(color)-.5],
                     plot_width=kwargs.get('plot_width', 200),
                     plot_height=kwargs.get('plot_height', 400),
                     title="Legend",
                     title_text_font_size=kwargs.get('title_text_font_size',
                                                     12),
                     x_axis_type=None, y_axis_type=None)
        lgd.square(x=1, y='y', size=20, color=color, source=source)
        lgd.text(x=1.1, y='y', text='text', source=source)
        lgd.grid.grid_line_color = None
        plist.append(lgd)
    if both:
        return gridplot([plist])
    else:
        return gridplot([plist])


def make_gridplot(y, df, x="i", text="None",
                  kwtext={'text_font_size': "6"}, groups=[],
                  ncol=1, total_plot_width=1200,
                  xaxis={'axis_label': "", 'major_label_orientation': np.pi/3},
                  yaxis={'axis_label': "", 'major_label_orientation': 1},
                  grid={'grid_line_color': "gray", 'grid_line_alpha': 0.2},
                  share_x_range=False, share_y_range=False,
                  tooltips=[], qc=None, abline={}, **kwargs):
    """Make a gridplot.

    Args:
      y (str): column name for y variable
      df (py:class:`~pandas.DataFrame`): data frame
      x (str): column name for x variable
      kwtext (dict): args passed to figure text object,
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
    if groups:
        grouped = df.groupby(groups)
    else:
        grouped = [("all", df)]
    first = True
    plist = []
    kwargs['plot_width'] = int(total_plot_width / ncol)
    kwargs['plot_height'] = int(total_plot_width / ncol)
    for name, data in grouped:
        source = ColumnDataSource(data)
        p = figure(title=name, **kwargs)
        if text:
            p.text(x=x, y=y, source=source, **kwtext)
        else:
            p.circle(x=x, y=y, source=source)
        if abline:
            x0 = abline.get('intercept', min(data[x]))
            x1 = max(data[x] * abline.get('pad', 1.1))
            y0 = abline.get('intercept', x0)
            y1 = ((x1-x0)*abline.get('slope', 1) + y0)
            p.line(x=[x0, x1], y=[y0, y1],
                   color=abline.get('color', 'red'))
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
