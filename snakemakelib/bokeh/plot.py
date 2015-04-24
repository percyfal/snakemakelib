# Copyright (C) 2015 by Per Unneberg
import os
import inspect
import numpy as np
from snakemakelib.config import sml_base_path
from bokeh.plotting import figure, Figure
from snakemakelib.log import LoggerManager

smllogger = LoggerManager().getLogger(__name__)

# Move elsewhere?
class QCArgs(object):
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
                xaxis = {'axis_label' : "", 'major_label_orientation' : np.pi/3},
                yaxis = {'axis_label' : "", 'major_label_orientation' : 1},
                grid = {'grid_line_color' : None, 'grid_line_alpha' : 1.0},
                tooltips=[], qc=None, **kwargs):
    """Make a scatter plot"""
    figmembers = [x[0] for x in inspect.getmembers(Figure)]
    circle_kwargs_keys = list(set(list(kwargs.keys())).difference(set(figmembers)))
    circle_kwargs = {k:kwargs.pop(k) for k in circle_kwargs_keys}
    # Get a reference to a figure
    fig = figure (**kwargs)
    # Add qc cutoff if required
    if not qc is None:
        getattr(fig, qc.plot_type)(x=qc.x, y=qc.y, **qc.kwargs)
    # Add scatter points - check if list of strings, if present in source we do several
    if set(y) <= set(circle_kwargs['source'].column_names):
        color = ["black"]
        if 'color' in circle_kwargs.keys():
            color = circle_kwargs.pop('color')
        for (yy, c) in zip(y, color):
            fig.circle(x=x, y=yy, color=c, **circle_kwargs)
    else:
        fig.circle(x=x, y=y, **circle_kwargs)
    # Apparently fig.xaxis and friends are lists so need to loop here. 
    for k in xaxis.keys():
        [setattr(x, k, xaxis[k]) for x in fig.xaxis]
    for k in yaxis.keys():
        [setattr(y, k, yaxis[k]) for y in fig.yaxis]
    for k in grid.keys():
        [setattr(x, k, grid[k]) for x in fig.grid]
    # NB: currently assume it is a dictionary
    for tt in tooltips:
        h = fig.select(dict(type=tt['type']))
        h.tooltips = tt['tips']
    return (fig)
    