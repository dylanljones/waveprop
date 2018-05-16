# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module contains 2-dimensional wrapper-classes of matplotlib-plots for easier plotting of data

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from .plot_utils import color_gaps


class Plot:
    """ Wrapper for matplotlib-plot for easy plotting

    Uses matplotlib for plotting package relevant data
    """

    E_LABEL = "$E \ [E_h]$"  # energy label
    V_LABEL = "$V \ [E_h]$"  # potential label
    X_LABEL = "$x \ [a_0]$"  # distance label

    default_font = {"size": 11}  # font for default plot
    double_font = {"size": 13}  # font for double plot
    quad_font = {"size": 15}  # font for quad plot

    default_width = 16  # default with of plot (cm)
    default_ratio = 3 / 5  # default ratio of plot
    default_linewidth = 1.5  # defauöt linewidth of plot
    default_dpi = 600  # default dpi for saving plot

    def __init__(self, title=None, subplots=None, height_ratios=None, width_ratios=None, sharex=None, sharey=None,
                 show=True, font=None, dpi=None, linewidth=None):
        """ Create and configure Plot instance

        Parameters
        ----------
        title : str
            title of plot
        subplots : array_like of int
            Number of rows/columns of the subplot grid.
        height_ratios : array_like of int
            ratios of cubplots in y-direction
        width_ratios : array_like of int
            ratios of cubplots in x-direction
        sharex : str
            Controls sharing of properties among x
        sharey : str
            Controls sharing of properties among y
        show : bool
            True if plot should be shown
        font : dict
            font to use in plot
        dpi : int
            dpi for saving plot
        linewidth : float
            default linewidth of plot
        """
        font = self.default_font if font is None else font
        self.set_font(font)
        self._dpi = None
        self._default_linewidth = linewidth if linewidth is not None else self.default_linewidth
        self._default_linestyle = "-"
        self._current_ax = 0

        if not show:
            plt.ioff()
        if subplots is None:
            fig, ax = plt.subplots()
            axs = [ax]
        else:
            height_ratios = [1] * subplots[0] if height_ratios is None else height_ratios
            width_ratios = [1] * subplots[1] if width_ratios is None else width_ratios
            sharex = "none" if sharex is None else sharex
            sharey = "none" if sharey is None else sharey
            fig, axs = plt.subplots(*subplots, sharex=sharex, sharey=sharey,
                                    gridspec_kw={"height_ratios": height_ratios, "width_ratios": width_ratios})
        self._fig, self._axs = fig, axs
        self.fig_setup(dpi=dpi)
        if title:
            self.set_title(title)

    @property
    def ax(self):
        """matplotlib.axes.Axes : currently set axes"""
        return self._axs[self._current_ax]

    @property
    def fig(self):
        """matplotlib.figure.Figure : Figure object"""
        return self._fig

    def fig_setup(self, width=None, ratio=None, dpi=None):
        """Set the figure appearence

        Parameters
        ----------
        width : float
            width of plot instance
        ratio : float
            ratio of plot instance
        dpi : int
            dpi for saving plot
        """
        width = width if width is not None else self.default_width
        ratio = ratio if ratio is not None else self.default_ratio
        dpi = dpi if dpi is not None else self.default_dpi

        w = width/2.54
        h = w*ratio
        self._fig.set_size_inches(w, h)
        self._dpi = dpi

    def set_default(self, width=1, ls="-"):
        """Set plot default settings

        Parameters
        ----------
        width : float
            default linewidth
        ls : str
            default linestyle
        """
        self._default_linewidth = width
        self._default_linestyle = ls

    def set_ax(self, i):
        """Set current axes

        Parameters
        ----------
        i : int
            index of axes
        """
        n = len(self._axs)
        if n == 1:
            self._current_ax = 0
        if i > n - 1:
            i = n - 1
        self._current_ax = i

    """ ======================================== Script plots ======================================================="""

    def script_layout(self, ratio=None, pad=None, spacing=None):
        """Set layout for latex-script

        Parameters
        ----------
        ratio : float
            ratio of plot instance
        pad : int
            padding to sides of plot instance
        spacing : int
            spacing between suplots
        """
        ratio = self.default_ratio if ratio is None else ratio
        pad = 2 if pad is None else pad
        spacing = 2 if spacing is None else spacing
        self.fig_setup(ratio=ratio)
        self.tight_layout(pad=pad, spacing=spacing)

    """ =========================================== Plotting ========================================================"""

    @classmethod
    def fast_plot(cls, x, y):
        """Class-Method for quick-plotting data

        Parameters
        ----------
        x : array_like
            x-data of plot
        y : array_like
            y-data of plot

        Returns
        -------
        plt : Plot
        """
        plot = cls()
        plot.plot(x, y)
        plot.show()
        return plot

    def plot(self, x, y, col=None, ls=None, width=None, label=None, **kwargs):
        """ Plot data

        Parameters
        ----------
        x : array_like
            x-data of plot
        y : array_like
            y-data of plot
        col : matplotlib color, optional
            color of line
        ls : str, optional
            linestyle
        width : float, optional
            linewidth
        label : str, optional
            label of data
        """
        width = self._default_linewidth if width is None else width
        ls = self._default_linestyle if ls is None else ls
        self.ax.plot(x, y, label=label, color=col, ls=ls, linewidth=width, **kwargs)

    def plot_curve(self, curve, col=None, label="", size="3", linewidth=None, plot_fit=True, fit_col=None, ls=None):
        """Plot data of curve object

        Parameters
        ----------
        curve : Curve
            curve object to be plotted
        col : matplotlib color
            color for plot
        label : str, optional
            label of data
        size : str, optional
            size of marker
        linewidth : float, optional
            width of line of fit
        plot_fit : bool, optional
            plot fit of curve if true
        fit_col : matplotlib color, optional
            color of fit-line
        ls : str, optional
            linestyle of fit-line
        """
        linewidth = self._default_linewidth if linewidth is None else linewidth
        if label == "":
            label = curve.label
        if curve.has_fit and plot_fit:
            if not fit_col:
                fit_col = col
            self.ax.plot(curve.fit_x, curve.fit_y, color=fit_col, linewidth=linewidth, ls=ls)
        self.plot_points(*curve.data, label, col, size)

    def plot_fit(self, curve, col=None, ls=None, width=None, label=None):
        """Plot fit of curve object

        Parameters
        ----------
        curve : Curve
            curve object of wich fit will be plotted
        col : matplotlib color, optional
            color for plot
        ls : str, optional
            linestyle
        width : float
            width of line of fit
        label : str, optional
            label of data
        """
        x, y = curve.fit_x, curve.fit_y
        self.plot(x, y, col=col, ls=ls, width=width, label=label)

    def plot_points(self, points, errs=None, label="", col=None, size="3"):
        """ Plot data points

        Parameters
        ----------
        points : array_like of array_like
            data of points
        errs : array_like of array_like, optional
            errors of points
        label : str, optional
            label of data
        col : matplotlib color, optional
            color for plot
        size : str, optional
            size of marker
        """
        x = np.asarray(points[0])
        y = np.asarray(points[1])
        x_err = np.asarray(errs[0])
        y_err = np.asarray(errs[1])

        x_and_y_errors = (x_err.shape[0] > 0) and (y_err.shape[0] > 0)
        only_x_errors = (x_err.shape[0] > 0) and (y_err.shape[0] == 0)
        only_y_erros = (x_err.shape[0] == 0) and (y_err.shape[0] > 0)

        if x_and_y_errors:
            self.ax.errorbar(x, y, x_err=x_err, yerr=y_err, fmt="o", markersize=size, label=label, capsize=2, color=col)
        elif only_x_errors:
            self.ax.errorbar(x, y, x_err=x_err, fmt="o", markersize=size, label=label, capsize=2, color=col)
        elif only_y_erros:
            self.ax.errorbar(x, y, yerr=y_err, fmt="o", markersize=size, label=label, capsize=2, color=col)
        else:
            self.ax.errorbar(x, y, fmt="o", markersize=size, label=label, capsize=2, color=col)

    def plot_values(self, points, line=0.5, marker="o", size="3", label="", col=None):
        x = np.asarray(points[0])
        y = np.asarray(points[1])
        self.ax.plot(x, y, marker=marker, markersize=size, label=label, color=col, linewidth=line)

    def plot_transmission(self, data, width=None, col=None, ls=None, label=None):
        """Quick plot transmission data

        Parameters
        ----------
        data : array_like of array_like
            transmission data
        width : float, optional
            linewidth of data
        col : matplotlib color, optional
            color of plot
        ls : str, optional
            linestyle
        label : str, optional
            label of data
        """
        x, y = data
        self.plot(x, y, label=label, col=col, ls=ls, width=width)
        self.set_labels(self.E_LABEL, "T")
        self.set_limits(xlim=[x[0], x[-1]], ylim=[0, 1])

    def show_band_strukture(self, model):
        """ Color band-gaps of given model on plot instance

        Parameters
        ----------
        model : Model
        """
        color_gaps(self, model)
        self.show()

    """ =========================================== Settings ========================================================"""

    @staticmethod
    def set_font(font):
        """ Set font of plot instance

        Parameters
        ----------
        font : dict
            font dict
        """
        matplotlib.rc('font', **font)

    def set_ticks_scale(self, axis, tick_scale):
        """ Set tick scale of given axis

        Parameters
        ----------
        axis : str
            name of axis
        tick_scale : float
            scale to set axis-ticking
        """
        if axis == "x":
            xlim = self.ax.get_xlim()
            ticks = np.arange(xlim[0], xlim[1]*1.1, tick_scale)
            self.ax.set_xticks(ticks)
        if axis == "y":
            ylim = self.ax.get_ylim()
            ticks = np.arange(ylim[0], ylim[1]*1.1, tick_scale)
            self.ax.set_yticks(ticks)

    def set_ticks_number(self, axis, n):
        """ Set number of ticks of given axis

        Parameters
        ----------
        axis : str
            name of axis
        n : int
            number of ticks
        """
        if axis == "x":
            xlim = self.ax.get_xlim()
            ticks = np.linspace(xlim[0], xlim[1]*1.1, n)
            self.ax.set_xticks(ticks)
        if axis == "y":
            ylim = self.ax.get_ylim()
            ticks = np.linspace(ylim[0], ylim[1]*1.1, n)
            self.ax.set_yticks(ticks)

    def config_axis(self, axis, ticks, labels=None):
        """ Set ticks and/or labels of chosen axis

        Parameters
        ----------
        axis : str
            name of axis
        ticks : array_like of float
            tick values
        labels : array_like of str, optional
            tick labels
        """
        if axis == "x":
            self.ax.set_xticks(ticks)
            if labels:
                self.ax.set_xticklabels(labels, minor=False)
        if axis == "y":
            self.ax.set_yticks(ticks)
            if labels:
                self.ax.set_yticklabels(labels, minor=False)

    def remove_spines(self, spines):
        """ Removes spines of plot for drawing

        Parameters
        ----------
        spines : array_like of str
            spines to remove
        """
        for spine in spines:
            self.ax.spines[spine].set_visible(False)

    def second_axis(self, axis, ticks, labels, same_scale=True):
        """ Add ticks on the other side of normal axis

        Parameters
        ----------
        axis : str
            name of axis
        ticks : array_like of float
            tick values
        labels : array_like of str
            tick labels
        same_scale : bool, optional
            use same scale on both axis

        Returns
        -------
        ax2 : matplotlib.axes.Axes
            second axis object
        """
        if axis == "y":
            ax2 = self.ax.twinx()
            if same_scale:
                ax2.set_ylim(*self.ax.get_ylim())
            ax2.set_yticks(ticks)
            ax2.set_yticklabels(labels, minor=False)
            return ax2
        if axis == "x":
            ax2 = self.ax.twiny()
            if same_scale:
                ax2.set_xlim(*self.ax.get_xlim())
            ax2.set_xticks(ticks)
            ax2.set_xticklabels(labels, minor=False)
            return ax2

    def set_title(self, title=""):
        """Sets title of plot instance

        Parameters
        ----------
        title : str
            title of plot
        """
        self.ax.set_title(title)

    def set_labels(self, xlabel=None, ylabel=None, xunit=None, yunit=None):
        """ Sets labels of plot instance

        Parameters
        ----------
        xlabel : str, optional
            label of x-axis
        ylabel : str, optional
            label of y-axis
        xunit : str, optional
            unit of x-axis
        yunit : str, optional
            unit of y-axis
        """
        if xlabel:
            label = xlabel + " [" + xunit + "]" if xunit else xlabel
            self.ax.set_xlabel(label)
        if ylabel:
            label = ylabel + " [" + yunit + "]" if yunit else ylabel
            self.ax.set_ylabel(label)

    def set_limits(self, xlim=None, ylim=None):
        """ Sets axis-limits of plot instance

        Parameters
        ----------
        xlim : array_like of float
            range o fx-axis
        ylim : array_like of float
            range of y-axis
        """
        if xlim:
            self.ax.set_xlim(xlim[0], xlim[1])
        if ylim:
            self.ax.set_ylim(ylim[0], ylim[1])

    """ ======================================== Drawing and Text ==================================================="""

    def add_lines(self, x=None, y=None, col="0.5", width=1, ls="-", label=None):
        """ Add straight lines to plot

        Parameters
        ----------
        x : array_like of float or float, optional
            x-values for lines
        y : array_like of float or float, optional
            y-values for lines
        col : matplotlib color, optional
            color of line
        width : float, optional
            linewidth
        ls : str, optional
            linestyle
        label : str, optional
            label of line
        """
        if x is not None:
            if not isinstance(x, list):
                x = [x]
            for i, _x in enumerate(x):
                if label and i == 0:
                    self.ax.axvline(_x, linestyle=ls, color=col, linewidth=width, label=label)
                else:
                    self.ax.axvline(_x, linestyle=ls, color=col, linewidth=width)
        if y is not None:
            if not isinstance(y, list):
                y = [y]
            for _y in y:
                self.ax.axhline(_y, linestyle=ls, color=col, linewidth=width)

    def fill(self, x, y, col="0.8", label=None):
        """ Fill range of plot with color

        Parameters
        ----------
        x : array_like of float
            start- and end-x-value for filling
        y : array_like of float
            start- and end-y-value for filling
        col : matplotlib color, optional
            color for fill
        label : str, optional
            label oof fill-region
        """
        self.ax.fill_between([x[0], x[1]], [y[1], y[1]], y[0], facecolor=col, alpha=0.5, label=label)

    def fill_between(self, x, y, const, col="0.8", label=None):
        """ Wrapper for matplotlib.Plot.fill_between

        Parameters
        ----------
        x : array_like
            x-data for fill
        y : array_like
            y-dta for fill
        const : float
            fill to this constant
        col : matplotlib color, optional
            color for fill
        label : str, optional
            label oof fill-region
        """
        self.ax.fill_between(x, y, const, facecolor=col, alpha=0.5, label=label)

    def text(self, pos, text, valign="center", halign="center"):
        """ Adds text to plot

        Parameters
        ----------
        pos : array_like of flota
            position of text
        text : str
            text to add
        valign : str, optional
            verticle alignment setting
        halign : str, optional
            horizontal alignment setting
        """
        x, y = pos
        self.ax.text(x, y, text, verticalalignment=valign, horizontalalignment=halign)

    def add_arrow_width(self, x0, x1, y, text, text_size=12, arrow_size=15, offset=0.5):
        """ Adds an horizontal arrow for describing a length

        Parameters
        ----------
        x0 : float
            start of arrow
        x1 : float
            end of arrow
        y : float
            y-value of arrow
        text : str
            text ontop of arrow
        text_size : int, optional
            size of text
        arrow_size : int, optional
            size of arrow
        offset : float, optional
            y-offset of text to arrow
        """
        x_c = x0 + (x1 - x0) / 2
        self.ax.text(x_c, y + offset, text, fontsize=text_size,
                     verticalalignment='center', horizontalalignment='center')
        self.ax.annotate("", xy=(x0, y), xytext=(x1, y), arrowprops=dict(arrowstyle='<->'), size=arrow_size)

    """ ============================================= Other ========================================================="""

    @staticmethod
    def close():
        """Wrapper-method for closing plot"""
        plt.close("all")

    def show(self, tight_layout=False):
        """Wrapper-method for showing plot

        Parameters
        ----------
        tight_layout : bool
            if true, set tight layout here for faster plotting
        """
        if tight_layout:
            self.tight_layout()
        plt.show()

    def save(self, path, dpi=None, box=None, file_type=".png"):
        """Wrapper method for saving plots

        Parameters
        ----------
        path : str
            location for saving
        dpi : int, optional
            dpi setting for saving
        box : str, optional
            bos setting for saving image
        file_type : str, optional
            filetype of image
        """
        name = path + file_type
        if dpi is None:
            dpi = self._dpi
        self._fig.savefig(name, dpi=dpi, bbox_inches=box)

    def tight_layout(self, pad=2, spacing=2):
        """ Wrapper method for tight-layout

        Parameters
        ----------
        pad : int, optional
            padding to border of plot
        spacing : int, optional
            spacing between subplots
        """
        self.fig.tight_layout(pad=pad, w_pad=spacing, h_pad=spacing)

    def legend(self, loc=None):
        """ Adds legend to plot

        Parameters
        ----------
        loc : str, optional
            location indicator for legend
        """
        if loc:
            location = loc
        else:
            location = 0
        self.ax.legend(loc=location)

    @staticmethod
    def color(index):
        """Wrapper-method for getting color by index"""
        return plt.get_cmap("tab10")(index)

    @staticmethod
    def linestyle(index):
        """Wrapper-method for getting linestyle by index"""
        styles = ["-", "--", "-.", ":"]
        return styles[index]


class ErrorPlot(Plot):
    """ Subclass Instance of plot to display data and error in subplot

    Data is displayed in top subplot, data in lower subplot
    """

    def __init__(self, title=None, height_ratio=3, font=None, dpi=600):
        """ Create subclass instance

        Parameters
        ----------
        title : str
            title of plot
        height_ratio : int
            ratio of main and error subplot
        font : dict
            font to use in plot
        dpi : int
            dpi for saving plot
        """
        super(ErrorPlot, self).__init__(title, font=font, subplots=(2, 1), height_ratios=[height_ratio, 1],
                                        sharex="all", dpi=dpi)
        self._error_limits = None
        self.set_ax(1)
        self.add_lines(y=0)
        self.set_ax(0)

    def set_labels(self, xlabel=None, ylabel=None, xunit=None, yunit=None):
        """ Sets labels of plot. Adds delta of error automaticly

        Parameters
        ----------
        xlabel : str, optional
            label of x-axis
        ylabel : str, optional
            label of y-axis
        xunit : str, optional
            unit of x-axis
        yunit : str, optional
            unit of y-axis
        """
        if xlabel:
            label = xlabel + " [" + xunit + "]" if xunit else xlabel
            self._axs[1].set_xlabel(label)
        if ylabel:
            label = ylabel + " [" + yunit + "]" if yunit else ylabel
            self._axs[0].set_ylabel(label)
            self._axs[1].set_ylabel("$\Delta$ " + label)

    def set_error_limit(self, ylim):
        """ Sets limits of error-subplot

        Parameters
        ----------
        ylim : array_like of float
            limits of error-plot
        """
        self.set_ax(1)
        self.set_limits(ylim=ylim)
        self.set_ax(0)

    def set_error_ticks(self, ticks, labels=None):
        """ Sets ticks and/or labels of error-subplot

        Parameters
        ----------
        ticks : list of float
            tick values
        labels : list of str
            tick labels
        """
        self.set_ax(1)
        self.config_axis("y", ticks=ticks, labels=labels)
        self.set_ax(0)

    def plot_curve(self, curve, col=None, label="", size="3", linewidth=None, plot_fit=True, fit_col=None, ls=None):
        """Wrapper for plot_curve-method of parent class"""
        self.set_ax(0)
        super(ErrorPlot, self).plot_curve(curve, col, label, size, linewidth, plot_fit)

    def plot_difference(self, curve, func, col=None, width=None, ls=None):
        """ Plots the difference of given data to given method in error-subplot

        Parameters
        ----------
        curve : Curve
            curve hoding data
        func : method
            method to calculate distance from data
        col : matplotlib color
            color of plot
        width : float
            linewidth
        ls : str
            linestyle
        """
        self.set_ax(1)
        x = curve.x
        delta = curve.y - np.vectorize(func)(x)
        self.plot(x, delta, col=col, width=width, ls=ls)

    def _plot_error(self, curve, col, offset=0.2):
        """ Plots the difference between fit data and data

        Parameters
        ----------
        curve : Curve
            curve object holding data and fit
        col : matplotlib color
            color of plot
        offset : float
            offset of plotting fit
        """
        offset += 1
        x = curve.x
        y_err = curve.y - curve.f(x)
        ymin, ymax = curve.error_limits

        self.set_ax(1)
        self._error_limits = offset * ymin, offset * ymax
        self.plot(x, y_err, col=col)
        self.set_ax(0)
