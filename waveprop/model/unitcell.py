# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module contains the Unitcell object

"""
import numpy as np


class Cell:
    """ Unit Cell of one-dimensional crystal modelled with rectangular potential barriers

    Each cell with size a is built up by a potential barrier with strength v and width d followed by
    some free space

      __      .        Cell n       .   Cell n+1   .
        |     . ___________         . _______      .
        |     .|           |        .|       |     . ____
        |     .|           | v      .|       |     .|
        |_____.|<--delta-->|________.|       |_____.|
              .<------- a --------->.              .
    """

    _cell_resolution = 100  #: Default cell resolution for plotting

    def __init__(self, v=10, a=1.0, d=0.8):
        """Create Cell Object

        Parameters
        ----------
        v : float
            Potential strength
        a : float
            Cell size
        d : float
            Barrier size
        """
        self._v = float(v)
        self._d = float(d)
        self._a = float(a)

    @property
    def params(self):
        """list of float: List of cell parameters v, a and d"""
        return [self._v, self._a, self._d]

    @property
    def a(self):
        """float: size of cell"""
        return round(self._a, 4)

    @property
    def v(self):
        """float: Potential strength of barrier"""
        return self._v

    @v.setter
    def v(self, v):
        self._v = float(v)

    @property
    def d(self):
        """float: Width of potential barrier"""
        return self._d

    @property
    def fr(self):
        """float: Length of free space"""
        return self._a - self._d

    def latex_str(self):
        """str:  Generates cell info in latex code"""
        return r"$V_s={:.1f}$, $a_s={:.1f}$, $d_s={:.1f}$".format(self.v, self.a, self.d).replace(".", ",")

    def build(self, x0, resolution=None):
        """Build potential data for plotting the cell

        Parameters
        ----------
        x0 : float
            Starting location of build
        resolution : int
            Number of points for building cell

        Returns
        -------
            potential_data : list of ndarray
        """
        if resolution is None:
            resolution = self._cell_resolution
        x_values = np.linspace(x0, x0 + self.a, resolution)
        v_values = []

        if self.fr == 0:
            v_values = [self.v]*len(x_values)

            v_values[-1] = 0
        else:
            for _x in x_values:
                x = _x - x0
                v = self.v if x <= self.d else 0
                v_values.append(v)
        v_values[0] = 0
        return [x_values, v_values]

    def __str__(self):
        out = "v={:.1f}, ".format(self._v)
        out += "a={:.2f}, ".format(self.a)
        out += "d={:.2f}".format(self._d)
        return out
