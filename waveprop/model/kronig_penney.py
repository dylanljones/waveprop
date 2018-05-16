# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module contains the Kronig-Penney-Model object

"""
import numpy as np
import cmath
from ..calculation import TransferMatrix, constants
from .unitcell import Cell
from .bands import BandFinder, bloch_vector, right_equation_from_e, get_gaps


class KronigPenney:
    """ Model of Kronig-Penney-Model

    Modelling of one-dimensional Kronig-Penney-Model consisting of an infinite amount
    of periodically aranged unitcells. Set Energy of System to calculate all attributes.
    """

    def __init__(self, cell=None, bands=3):
        """Create Model

        Parameters
        ----------
        cell : Cell
            Unitcell of model
        bands : int
            Number of bands to calculate
        """
        if cell is None:
            cell = Cell()

        self._cell = cell
        self._transferMatrix = TransferMatrix()
        self._bands = self._get_bands(bands)
        self._e = 0
        self._set()

    @classmethod
    def from_parameters(cls, v, a, d):
        """Create Object without Cell Object

        Parameters
        ----------
        v : float
            Potential strength
        a : float
            Cell size
        d : float
            Barrier size

        Returns
        -------
        Object: KronigPenney
        """
        cell = Cell(v, a, d)
        return cls(cell)

    def set_energy(self, e):
        """Set energy of particle moving through Model

        Parameters
        ----------
        e : float
            Energy of particle
        """
        q = np.sqrt(2 * constants.m * e) / constants.hbar
        kappa = cmath.sqrt(2 * constants.m * (self._cell.v - e)) / constants.hbar
        k = bloch_vector(e, self._cell)
        eq_r = right_equation_from_e(e, self._cell)
        self._e = e
        self._set(q, kappa, k, eq_r)

    @property
    def cell(self):
        """Cell: Unitcell used in model"""
        return self._cell

    @property
    def q(self):
        """float: Wavevector in free space"""
        return self._q

    @property
    def kappa(self):
        """float: Wavevector in barrier"""
        return self._kappa

    @property
    def k(self):
        """float: Bloch vector of model"""
        return self._k

    @property
    def M(self):
        """ndarray: Transfermatrix of the model"""
        return self._transferMatrix.M

    @property
    def t(self):
        """float: Transmission koefficient of model"""
        return self._transferMatrix.t

    @property
    def bands(self):
        """list of float : Energy bands of model"""
        return self._bands

    @property
    def gaps(self):
        """list of float : Energy gaps of model"""
        return get_gaps(self.bands)

    @property
    def bloch_equation(self):
        """float: Right side of transcendent equation"""
        return self._right_equation

    def _get_bands(self, n):
        """ Calculate the first n energy bands of the model using the BandFinder object

        Parameters
        ----------
        n : int
            Number of bands

        Returns
        -------
        bands list of float
            Energy bands of model
        """
        bf = BandFinder(self._cell, n)
        return bf.bands

    def _set(self, q=0, kappa=0, k=0, eq_r=0):
        """Sets all values of the object

        Parameters
        ----------
        q : float
            Wavevector in free space
        kappa : float
            Wavevector in barrier
        k : float
            Bloch vector of model
        eq_r : float
            Right side of transcendent equation
        """
        self._q = q
        self._kappa = kappa
        self._k = k
        self._right_equation = eq_r
        if np.imag(self._k) == 0:
            self._transferMatrix.single(q, kappa, self._cell)
        else:
            self._transferMatrix.null()

    def build(self, x0, n, cell_resolution=None):
        """Build potential data for plotting the model

        Parameters
        ----------
        x0 : float
            Starting location of build
        n :
            int
            Number of cells in build
        cell_resolution : int
            Number of points in each cell

        Returns
        -------
            potential_data : list of ndarray
        """
        x_values, v_values = np.asarray([]), np.asarray([])
        for _ in range(n):
            _x, _v = self.cell.build(x0, cell_resolution)
            x_values = np.append(x_values, _x)
            v_values = np.append(v_values, _v)
            x0 += self.cell.a
        return [x_values, v_values]

    def __str__(self):
        out = str(self.cell) + "\n"
        return out

    """
    def state_density(self):
        e = self._e
        n = 200
        e_values = np.linspace(e*0.99, e*1.01, n)
        k_values = []
        for e in e_values:
            k = np.real(bloch_vector(e, self.cell))
            k_values.append(k)
        grad = np.gradient(k_values)[int(n/2)-1]
        return 1/(2*np.pi)*abs(grad)
    """