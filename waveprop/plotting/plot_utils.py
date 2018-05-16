# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module contains helper-methods for plotting data

"""
import numpy as np


def color_gaps(plot, model, ylim=None, label=False, color_sample_gaps=True):
    """ Colors the energies of the band-gaps of the model

    Parameters
    ----------
    plot : matplotlib object
        plot to add bandgaps
    model : Model
        model to get band-data from
    ylim : array_like, optional
        y-range to color gaps,
        default is 0-1
    label : bool, optional
        if true, label the gap-energies,
        default: False
    color_sample_gaps: bool, optional
        if true, color gap-energies of ordered sample,
        default: True
    """
    col = 0.8
    if ylim is None:
        ylim = [0, 1]
    for i, g in enumerate(model.gaps):
        if label and i == 0:
            plot.fill(g, ylim, col=str(col), label="Leitungen")
        else:
            plot.fill(g, ylim, col=str(col))
    sample_gaps = model.sample_gaps
    if color_sample_gaps and sample_gaps:
        extra_gaps = []
        for band in model.bands:
            for sample_gap in sample_gaps:
                if any([model.in_band(x, band) for x in sample_gap]):
                    start = sample_gap[0]
                    if start < band[0]:
                        start = band[0]
                    end = sample_gap[1]
                    if end > band[1]:
                        end = band[1]
                    extra_gaps.append([start, end])
        for i, g in enumerate(extra_gaps):
            if label and i == 0:
                plot.fill(g, ylim, col=(1, col-0.1, col-0.1, 0.2), label="Streubereich")
            else:
                plot.fill(g, ylim, col=(1, col-0.1, col-0.1, 0.2))


def plot_model(model, font=None, cell_res=None, n_cells=None):
    """ Plot the potential data of model

    Parameters
    ----------
    model : Model
        model to plot
    font : dict, optional
        font to use in plot
    cell_res : int, optional
        number of points in one cell
    n_cells : int
        number of cells to show for each lead

    Returns
    -------
    plt : matplotlib object
    """
    from .plot_objects import Plot
    plt = Plot(font=font)
    ax = plt.ax
    if model.sample:
        plt.add_lines(x=[0, model.sample_length], col="0.7")

    xlim, ylim = add_model_plot(ax, model, cell_res, n_cells, True)
    plt.set_limits(xlim, ylim)

    if model.sample:
        plt.legend()
    plt.set_labels(plt.X_LABEL, plt.V_LABEL)
    return plt


def add_model_plot(ax, model, res, n_cells=None, label=True, sample_color=None):
    """ Add model-plot to matplotlib-figure

    Parameters
    ----------
    ax : matplotlib axis
        axis to add model-plot
    model : Model
        model to plot
    res : int
        number of points in one cell
    n_cells : int
        number of cells to show for each lead
    label : bool
        if True, label the Leads of the model
    sample_color : matplotlib color
        color to display the sample

    Returns
    -------
    lim_data : tuple of tuple
        x-limit and y-limit data of plot
    """
    alpha = 0.2
    lead_l, lead_r, sample = model.build(res, n_cells)
    _ymax = max(np.append(lead_l[1], sample[1]))
    ymax = _ymax + 0.15 * _ymax
    xmin, xmax = min(lead_l[0]), max(lead_r[0])

    if model.sample:
        color_lead = ymax * 2
    else:
        color_lead = None
    fill_color = "0.4"

    label = "Leads" if label else None
    _add_part(ax, lead_l, "black", fill_color, alpha, color_lead, label)
    _add_part(ax, lead_r, "black", fill_color, alpha, color_lead)
    _add_part(ax, sample, None, fill_color=sample_color, alpha=alpha)
    return (xmin, xmax), (0, ymax)


def _add_part(ax, data, color, fill_color=None, alpha=0.2, color_lead=None, label=None):
    """ Helper method to add parts of the model to the plot

    Parameters
    ----------
    ax : matplotlib axis
        axis to add model-plot
    data : array_like
        potential data of the part
    color : matplotlib color
        line color of potential
    fill_color : matplotlib color, optional
        color to fill the potential barriers
    alpha : float, optional
        alpha value for the filling color
    color_lead : bool, optional
        if True, mark the lead regions
    label : str, optional
        label of the part
    """
    if not fill_color:
        fill_color = color
    x, y = data[0], data[1]
    ax.plot(x, y, color=color)
    ax.fill_between(x, 0, y, facecolor=fill_color, alpha=alpha)
    if color_lead:
        ax.fill_between(x, color_lead, y, facecolor="0.8", alpha=0.5, label=label)
