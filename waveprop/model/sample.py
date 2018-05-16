# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module contains scattering sample objects:
    -Sample:            general Sample consisting of different unitcells
    -OrderedSample:     sample of n identical unicells
    -DisorderedSample:  randomly distributet sample

"""
import numpy as np
import random
from ..calculation import constants, TransferMatrix
from .unitcell import Cell
from .kronig_penney import BandFinder, get_gaps, bloch_vector


class Sample:
    """ General sample object

    general sample object consisting of multiple unitcells. other sample objects
    inherent from this
    """

    def __init__(self, cells):
        """Create sample object

        Parameters
        ----------
        cells : list of Cell
            list of unitcells in sample
        """
        self._cells = list(cells)
        self._transferMatrix = TransferMatrix()
        self._e = 0

    def set_energy(self, e):
        """Set energy of particle moving through Model

        Parameters
        ----------
        e : float
            Energy of particle
        """
        q = np.sqrt(2 * constants.m * e) / constants.hbar
        self._transferMatrix.multiple(q, e, self._cells)
        self._e = e

    @property
    def M(self):
        """ndarray: Transfermatrix of the model"""
        return self._transferMatrix.M

    @property
    def t(self):
        """float: Transmission koefficient of model"""
        return self._transferMatrix.t

    @property
    def cells(self):
        """list of Cell: list of cells in sample region"""
        return self._cells

    @cells.setter
    def cells(self, cells):
        self._cells = cells

    @property
    def n(self):
        """int: number of unitcells in sample"""
        return len(self._cells)

    @property
    def length(self):
        """float: total length of sample"""
        length = 0
        for cell in self._cells:
            length += cell.a
        return length

    def build(self, x0=0, cell_resolution=None):
        """Build potential data for plotting the sample

        Parameters
        ----------
        x0 : float
            Starting location of build
        cell_resolution : int
            Number of points in each cell

        Returns
        -------
            potential_data : list of ndarray
        """
        x_values, v_values = np.asarray([]), np.asarray([])
        for cell in self.cells:
            x, v = cell.build(x0, cell_resolution)
            x_values = np.append(x_values, x)
            v_values = np.append(v_values, v)
            x0 += cell.a
        return x_values, v_values

    def type(self):
        """type: get type of sample"""
        return type(self)

    def _to_string(self):
        """str: get info string of cell"""
        out = "N = {:d}, Length = {:.1f}".format(self.n, self.length)
        return out

    def __str__(self):
        out = "Sample: \n"
        return out


class OrderedSample(Sample):
    """ Ordered sample object

    sample consisting of n identical unitcells. Has band structure and bloch vector
    """

    def __init__(self, cell, n, n_bands=2):
        """Create ordered Sample

        Parameters
        ----------
        cell : Cell
            unitcell of sample
        n : int
            number of identical unitcells in sample
        n_bands : int
            number of sample-bands to calculate
        """
        self._cell_prototype = cell
        cells = [cell] * n
        super(OrderedSample, self).__init__(cells)
        self._bands = self.get_bands(n_bands)

    def new_sample(self, n, cell=None):
        """Generate new sample

        Parameters
        ----------
        n : int
            number of identical unitcells in sample
        cell : Cell
            unitcell of sample
        """
        if cell is None:
            cell = self._cell_prototype
        self.cells = [cell] * n

    @property
    def cell(self):
        """Cell: unitcell used in sample"""
        return self._cells[0]

    @property
    def bands(self):
        """list of float : Energy bands of sample"""
        return list(self._bands)

    @property
    def gaps(self):
        """list of float : Energy gaps of sample"""
        return get_gaps(self.bands)

    @property
    def k(self):
        """float: Bloch vector of sample"""
        return bloch_vector(self._e, self._cell_prototype)

    def get_bands(self, n=3):
        """ Calculate the first n energy bands of the model using the BandFinder object

        Parameters
        ----------
        n : int
            Number of bands

        Returns
        -------
        bands : list of list
            Energy bands of sample
        """
        bf = BandFinder(self._cell_prototype, n)
        return bf.bands

    def get_band_size(self, i):
        """Calculate size of i-th energy band of sample

        Parameters
        ----------
        i : int
            number of band

        Returns
        -------
        size : float
        """
        band = list(self.bands[i-1])
        return abs(band[1]-band[0])

    def get_gap(self, i):
        """list of float: i-th energy gap """
        return self.gaps[i]

    def get_gap_size(self, i):
        """Calculate size of i-th energy gap of sample

        Parameters
        ----------
        i : int
            number of gap

        Returns
        -------
        size : float
        """
        gap = self.gaps[i]
        return abs(gap[1]-gap[0])

    def __str__(self):
        out = "Ordered Sample: V = {}, ".format(self.cell.v)
        return out + self._to_string()


class DisorderedSample(Sample):
    """ Disordered sample object

    sample consisting of n random unitcells. The potential strength has a rectangular distribution
    with width w. This describes the disorder of the ssample.
    """

    def __init__(self, cell, w, n):
        """Crerate Disordered Sample

        Parameters
        ----------
        cell : Cell
            unitcell to base probability distribiution
        w : float
            strngth of disorder
        n : int
            number of unitcells in sample
        """
        self._cell_prototype = cell
        self._w = w
        cells = self._get_cells(self._w, n, cell)
        super(DisorderedSample, self).__init__(cells)

    @property
    def w(self):
        """float: strength of disorder in fraction of v of leads"""
        return self._w

    def new_sample(self, w, n, cell=None):
        """Generate new sample

        Parameters
        ----------
        w : float
            strngth of disorder
        n : int
            number of unitcells in sample
        cell : Cell
            unitcell to base probability distribiution
        """
        self._w = w
        if cell is None:
            cell = self._cell_prototype
        cells = self._get_cells(self._w, n, cell)
        self.cells = cells

    @staticmethod
    def _get_cells(w, n, cell):
        """Generate n unitcells with random barrier strengths v_i around v

        Parameters
        ----------
        w : float
            strngth of disorder
        n : int
            number of unitcells in sample
        cell : Cell
            unitcell to base probability distribiution

        Returns
        -------
        cells : list of Cell
            list of generated unitcell-objects
        """
        v0, a, d = cell.params
        delta = v0 * w
        cells = []
        for i in range(n):
            v = random.uniform(v0 - delta, v0 + delta)
            cell = Cell(v, a, d)
            cells.append(cell)
        return cells

    def __str__(self):
        out = "Disordered Sample: W = {}, ".format(self.w)
        return out + self._to_string()
