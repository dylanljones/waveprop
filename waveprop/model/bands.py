# -*- coding: utf-8 -*-
"""
Created on 8 May 2018
@author: Dylan Jones

This module contains mothods for calculating the band structure and the BandFinder object

"""

from cmath import acos, cos, cosh, sin, sinh, sqrt
from scipy.signal import argrelextrema
import numpy as np


def get_gaps(bands):
    """Calculate gaps from bands

    Parameters
    ----------
    bands : array_like of array_like of float
        Energy bands of model

    Returns
    -------
    gaps list of float
        Energy gaps of model
    """
    if not len(bands):
        return None
    b = bands[0]
    gaps = [[0, b[0]]]
    for i in range(0, len(bands) - 1):
        gaps.append([bands[i][1], bands[i + 1][0]])
    return gaps


def right_equation(q, kappa, cell):
    """Calculates the right side of the transzendent bloch-equation from wavevectors

    Parameters
    ----------
    q : float
        Wavevector in free space
    kappa : float
        Wavevector in barrier
    cell : Cell
        Unitcell of kronig-penney-model

    Returns
    -------
    eq_r : float
        Right side of transcendent equation
    """
    x = kappa * cell.d
    y = q * round(cell.fr, 2)
    return cosh(x)*cos(y) + (kappa**2 - q**2)/(2*kappa*q) * sinh(x)*sin(y)


def right_equation_from_e(e, cell):
    """Calculates the right side of the transzendent bloch-equation from energy

    Parameters
    ----------
    e : float
        Energy of the particle
    cell : Cell
        Unitcell of kronig-penney-model

    Returns
    -------
    eq_r : float
        Right side of transcendent equation
    """
    if e == 0 or e == cell.v:
        return np.nan
    q = sqrt(2 * e)
    kappa = sqrt(2 * (cell.v - e))
    return right_equation(q, kappa, cell)


def bloch_vector(e, cell):
    """Calculates the bloch vector of the system for the energy e

    Parameters
    ----------
    e : float
        Energy of the particle
    cell : Cell
        Unitcell of kronig-penney-model

    Returns
    -------
    k : float
        Bloch vector
    """
    eq_r = right_equation_from_e(e, cell)

    if abs(eq_r) < 1:
        k = np.arccos(eq_r) / cell.a
    elif abs(eq_r) == 1:
        k = 0
    else:
        k = acos(eq_r) / cell.a
    return k


def right_equation_limit(e, cell):
    """Calculates the energies where the right side of the transzendent bloch-equation equals +- 1

    Parameters
    ----------
    e : float
        Energy of the particle
    cell : Cell
        Unitcell of kronig-penney-model

    Returns
    -------
    limit: method
        Limit method of transzendent bloch-equation
    """
    d = np.real(right_equation_from_e(e, cell))
    limit = np.abs(np.abs(d) - 1)
    return limit


class BandFinder:
    """ Object to calculate band structure of kronig-penney-model

    Uses the right side of the transcendent bloch-equation f(E) to check
    if value is in the right range of
        |f(E)| <= 1.
    to get band edges, the energy values of |f(E)| = 1 are calculated iterative

    Class Attributes
    ----------
    bands : array_like of tuple of float
        band data containing lower and upper limit of bands

    """

    def __init__(self, cell, n_bands=2, n_x=1E4):
        """Create the object

        Parameters
        ----------
        cell : Cell
            Unitcell of kronig-penney-model
        n_bands : int
            Number of bands to calculate
        n_x : int
            Number of energy points per range
        """
        self._cell = cell
        self._n_x = int(n_x)
        self._f = np.vectorize(right_equation_limit)

        self._x = np.asarray([])
        self._limits = np.asarray([])
        self.bands = []

        self._steps = 5, 10, 100, 1000
        self._chunks = 0, 30, 100, 1000, 1E10

        self._calculate(n_bands)

    def lower_limit(self, index):
        """float: get lower limit of band"""
        return self.bands[index-1][0]

    def upper_limit(self, index):
        """float: get upper limit of band"""
        return self.bands[index-1][1]

    def _calculate(self, n_bands):
        """ Calculate band strukture and set class attributes

        Parameters
        ----------
        n_bands : int
            Number of bands to calculate
        """
        d = self._cell.d
        a = self._cell.a
        bands = []
        c = np.pi ** 2 / 2

        if d == 0:
            for i in range(n_bands):
                e0 = c * (2 * i / a) ** 2
                e1 = c * ((2 * i + 1) / a) ** 2
                band = [e0, e1]
                bands.append(band)

        elif d == a:
            v = self._cell.v
            for i in range(n_bands):
                e0 = c * (2 * i / a) ** 2 + v
                e1 = c * ((2 * i + 1) / a) ** 2 + v
                band = [e0, e1]
                bands.append(band)

        else:
            bands = self._search(n_bands)
        self.bands = bands

    def _search(self, n_bands):
        """ Search for band limits

        Parameters
        ----------
        n_bands : int
            Number of bands to calculate

        Returns
        -------
        bands list of float
            Energy bands of model
        """
        x1 = 0
        while len(self._limits) < (n_bands * 2):
            x0 = x1
            s = self._set_step_size(x0)
            x1 = x0 + s
            self._x = np.linspace(x0, x1, self._n_x)
            limits = self._get_limits()
            self._limits = np.append(self._limits, limits)
        return self._get_bands(n_bands)

    def _set_step_size(self, x0):
        """Set energy step size in dependence of energy range

        Parameters
        ----------
        x0 : float
            start of energy range

        Returns
        -------
        steps_size : float
            energy step size for calculation

        """
        if x0 == 0:
            return self._steps[0]
        for i in range(len(self._steps)):
            xmin, xmax = self._chunks[i], self._chunks[i + 1]
            bigger_as_xmin = xmin <= x0
            smaller_as_xmax = x0 < xmax
            if bigger_as_xmin and smaller_as_xmax:
                return self._steps[i]

    def _get_bands(self, n):
        """Build bands from limits of transcendent equation

        Parameters
        ----------
        n : int
            number of bands

        Returns
        -------
        bands list of float
            Energy bands of model
        """
        bands = []
        for i in range(0, n * 2, 2):
            bands.append([self._limits[i], self._limits[i + 1]])
        return bands

    def _get_limits(self):
        """ Find limits of transcendent equation f, where |f(e)| = 1

        Returns
        -------
        limits list of float
            energy levels for limit of transcendent equation
        """
        y = self._f(self._x, self._cell)
        y = self._clean_data(y)
        indexes = argrelextrema(y, np.less)
        mins = []
        for i in indexes:
            mins.append(self._x[i])
        return mins[0]

    def _clean_data(self, y_values):
        """Remove np.nan's from data arrays

        Parameters
        ----------
        y_values : list of float
            data to clean

        Returns
        -------
        cleaned_data : list of float
        """
        n = self._x.shape[0]
        x_out, y_out = [], []
        for i in range(n):
            x, y = self._x[i], y_values[i]
            if np.isnan(y):
                continue
            x_out.append(x)
            y_out.append(y)
        self._x = np.asarray(x_out)
        return np.asarray(y_out)
