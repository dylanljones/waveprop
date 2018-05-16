# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module contains a 3-dimensional wrapper-class of matplotlib-plots for easier plotting of data

"""
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D


class Plot3D:
    """ Wrapper for matplotlib-plot for easy 3D-plotting

    Uses matplotlib for plotting package relevant data
    """

    default_font = {"size": 11}  # font for default plot
    default_linewidth = 1.5  # # default linewidth for plotting
    default_ratio = 3 / 5  # default ratio of plot

    def __init__(self, font=None, dpi=600):
        """ Create 3d plot instance

        Parameters
        ----------
        font : dict
            font to use in plot
        dpi : int
            dpi for saving plot
        """
        font = self.default_font if font is None else font
        self.set_font(font)
        self.dpi = dpi

        fig = plt.figure()
        self._ax = Axes3D(fig)
        self._fig = fig
        self.fig_setup()

    def fig_setup(self, width=16, ratio=1, dpi=600):
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
        w = width/2.54
        h = w*ratio
        self._fig.set_size_inches(w, h)
        self.dpi = dpi

    @staticmethod
    def set_font(font):
        """ Set font of plot instance

        Parameters
        ----------
        font : dict
            font dict
        """
        matplotlib.rc('font', **font)

    @property
    def ax(self):
        """matplotlib.axes.Axes : axes object of plot"""
        return self._ax

    @property
    def fig(self):
        """matplotlib.figure.Figure : Figure object"""
        return self._fig

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

    def surface(self, x, y, z, col=None, alpha=None, cmap=None, cmap_limits=None):
        """ Plots 3D data as surface

        Parameters
        ----------
        x : array_like
        y : array_like
        z : array_like
        col : matplotlib color
            color of sufrace
        alpha : float
            alpha value of surface
        cmap : str
            color-map for surface
        cmap_limits : array_like of float
            limits for color-map
        """
        cmap_limits = [None, None] if cmap_limits is None else cmap_limits
        self.ax.plot_surface(x, y, z, lw=0.5, color=col, rstride=1, cstride=1, alpha=alpha,
                             cmap=cmap, vmin=cmap_limits[0], vmax=cmap_limits[1], antialiased=False)

    def wireframe(self, x, y, z, col=None, rstride=5, cstride=5, alpha=0.5):
        """ Plots 3D data as wireframe

        Parameters
        ----------
        x : array_like
        y : array_like
        z : array_like
        col : matplotlib color
            color of wireframe
        rstride : int
            number of wires in y-direction (rows)
        cstride : int
            number of wires in x-direction (collumns)
        alpha : float
            alpha value of wirefram
        """
        self.ax.plot_wireframe(x, y, z, color=col, rstride=rstride, cstride=cstride, alpha=alpha)

    def line(self, x, y, z, col=None, alpha=1):
        """ Plots 3D data as line

        Parameters
        ----------
        x : array_like
        y : array_like
        z : array_like
        col : matplotlib color
            color of line
        alpha : float
            alpha value of line
        """
        self.ax.plot(x, y, z, color=col, alpha=alpha)

    def contour_3d(self, x, y, z, n, col=None, alpha=None, cmap=None, cmap_limits=None):
        """ Plots 3D data as contours

        Parameters
        ----------
        x : array_like
        y : array_like
        z : array_like
        n : int
            number of contour-lines
        col : matplotlib color
            color of line
        alpha : float
            alpha value of line
        cmap : str
            color-map for contour
        cmap_limits : array_like of float
            limits for color-map
        """
        cmap_limits = [None, None] if cmap_limits is None else cmap_limits
        self.ax.contour3D(x, y, z, n, lw=0.5, color=col, rstride=1, cstride=1, alpha=alpha,
                          cmap=cmap, vmin=cmap_limits[0], vmax=cmap_limits[1], antialiased=False)

    """ =========================================== Settings ========================================================"""

    def set_labels(self, xlabel=None, ylabel=None, zlabel=None):
        """ Sets labels of plot instance

        Parameters
        ----------
        xlabel : str, optional
            label of x-axis
        ylabel : str, optional
            label of y-axis
        zlabel : str, optional
            label of z-axis
        """
        if xlabel:
            self.ax.set_xlabel(xlabel)
        if ylabel:
            self.ax.set_ylabel(ylabel)
        if zlabel:
            self.ax.set_zlabel(zlabel)

    def set_limits(self, xlim=None, ylim=None, zlim=None):
        """ Sets axis-limits of plot instance

        Parameters
        ----------
        xlim : array_like
            range o fx-axis
        ylim : array_like
            range of y-axis
        zlim : array_like
            range of z-axis
        """
        if xlim:
            self.ax.set_xlim(xlim[0], xlim[1])
        if ylim:
            self.ax.set_ylim(ylim[0], ylim[1])
        if zlim:
            self.ax.set_zlim(zlim[0], zlim[1])

    def config_axis(self, axis, ticks, labels=None):
        """ Set ticks and/or labels of chosen axis

        Parameters
        ----------
        axis : str
            name of axis
        ticks : array_like
            tick values
        labels : array_like of str
            tick labels
        """
        if axis == "x":
            self.ax.set_xticks(ticks)
            if labels:
                self.ax.set_xticklabels(labels, minor=False)
        elif axis == "y":
            self.ax.set_yticks(ticks)
            if labels:
                self.ax.set_yticklabels(labels, minor=False)
        elif axis == "z":
            self.ax.set_zticks(ticks)
            if labels:
                self.ax.set_zticklabels(labels, minor=False)

    """ ============================================= Other ========================================================="""

    def get_view(self):
        elev = self.ax.elev
        azim = self.ax.azim
        dist = self.ax.dist
        return azim, elev, dist

    def set_view(self, azim, elev=None, dist=None):
        self.ax.view_init(azim=azim, elev=elev)
        if dist:
            self.ax.dist = dist

    @staticmethod
    def close():
        """Wrapper-method for closing plot"""
        plt.close("all")

    @staticmethod
    def show():
        """Wrapper-method for showing plot"""
        plt.show()

    def save(self, path, dpi=None, box=None, file_type=".png"):
        """Wrapper method for saving plots

        Parameters
        ----------
        path : str
            location for saving
        dpi : int
            dpi setting for saving
        box : str
            bos setting for saving image
        file_type : str
            filetype of image
        """
        name = path + file_type
        if dpi is None:
            dpi = self.dpi
        self._fig.savefig(name, dpi=dpi, bbox_inches=box)

    def tight_layout(self, pad=2, spacing=2):
        """ Wrapper method for tight-layout

        Parameters
        ----------
        pad : int
            padding to border of plot
        spacing : int
            spacing between subplots
        """
        self.fig.tight_layout(pad=pad, w_pad=spacing, h_pad=spacing)

    def legend(self, loc=None):
        """ Adds legend to plot

        Parameters
        ----------
        loc : str
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
    def cmap(index):
        """Method for getting colormaps from index"""
        maps = ["coolwarm", "Blues"]
        return maps[index]
