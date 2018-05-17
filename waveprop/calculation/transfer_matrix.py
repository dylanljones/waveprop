# -*- coding: utf-8 -*-
"""
Created on 6 May 2018
@author: Dylan Jones

This module inplements the transfermatrix-method for calculating the transmission
through a scattering region

"""

import numpy as np
from numpy import exp, sinh, cosh, linalg
import cmath
from .constants import constants
from ..model import Cell


class TransferMatrix:
    """ transfermatrix object

    object for caclulating the transfermatrix and the transmission-koefficient
    of a scattering region
    """

    def __init__(self):
        """Create object and initialize transfermatrix"""
        self._M = np.zeros((2, 2))

    @property
    def M(self):
        """ndarray: transfermatrix of scattering region for the set energy"""
        return self._M

    @property
    def t(self):
        """float: transmission-koefficinet of scattering region for set energy"""
        m22 = self._M[1, 1]
        if m22 == 0:
            return 0
        else:
            return 1 / abs(m22) ** 2

    def null(self):
        """Set to an empty transfermatrix"""
        self._M = np.zeros((2, 2))

    def barrier(self, e, v, d):
        """Configure the transfermatrix of a rectangular potential barrier

        Parameters
        ----------
        e : float
            Energy of the particle
        v : float
            Potential strength of the barrier
        d : float
            Width of the barrier
        """
        v_in = cmath.sqrt(2 * constants.m * (v - e)) / constants.hbar
        v_out = cmath.sqrt(2 * constants.m * e) / constants.hbar
        self._M = self._transfer_matrix_barrier(v_out, v_in, d)

    def cells(self, e, cells):
        """Configure the transfermatrix of multiple unitcells for energy e

        Parameters
        ----------
        e : float
            Energy of the particle
        cells : array_like of Cell or Cell
            unitcell(s) of the scattering region
        """
        if isinstance(cells, Cell):
            cells = [cells]
        v_out = cmath.sqrt(2 * constants.m * e) / constants.hbar
        self.multiple(v_out, e, cells)

    def single(self, v_out, v_in, cell):
        """Configure the transfermatrix of a single unitcell

        Parameters
        ----------
        v_out : float
            Wavevector in free space
        v_in : float
            Wavevector in the potential barrier
        cell : Cell
            unitcell of scattering region
        """
        self._M = self._transfer_matrix_unitcell(v_out, v_in, cell)

    def multiple(self, v_out, e, cells):
        """Configure the transfermatrix of multiple unitcells

        Parameters
        ----------
        v_out : float
            Wavevector in free space
        e : float
            Energy of the particle
        cells : array_like of Cell
            unitcells of the scattering region
        """
        m_total = 1
        for cell in cells:
            v_in = cmath.sqrt(2 * constants.m * (cell.v - e)) / constants.hbar
            m = self._transfer_matrix_unitcell(v_out, v_in, cell)
            m_total = np.dot(m, m_total)
        self._M = m_total

    def diagonalize(self, base_matrix):
        """Represent the current transfermatrix in the base of the given matrix

        Parameters
        ----------
        base_matrix : ndarray
            Base matrix
        """
        q = linalg.eig(base_matrix)[1]
        m_diag = np.dot(linalg.inv(q), np.dot(self._M, q))
        self._M = m_diag

    def transmission_curve(self, xlim, cells, steps=1000):
        """Calculate transmission values for energies in a given range

        Parameters
        ----------
        xlim : array_like of float
            energy range for calculating transmission curve, consisting of
            the start and end value.
        cells : array_like of Cell
            unitcells of the scattering region
        steps : int
            number of energy levels to calculate

        Returns
        -------
        data : array_like of ndarray
            e and t data of the transmission curve
        """
        e_values = np.linspace(*xlim, steps)
        t_values = []
        for e in e_values:
            self.cells(e, cells)
            t_values.append(self.t)
        return e_values, t_values

    def _transfer_matrix_unitcell(self, v_out, v_in, cell):
        """Calculate the transfer matrix for a unitcell

        Parameters
        ----------
        v_out : float
            Wavevector in free space
        v_in : float
            Wavevector in the potential barrier
        cell : Cell
            unitcell

        Returns
        -------
        M : ndarray
        """
        if cell.fr == 0:
            m = self._transfer_matrix_barrier(v_out, v_in, cell.d)
        else:
            m_b = self._transfer_matrix_barrier(v_out, v_in, cell.d)
            m_f = self._transfer_matrix_free(v_out, cell.fr)
            m = np.dot(m_f, m_b)
        return m

    @staticmethod
    def _transfer_matrix_barrier(v_out, v_in, d):
        """Calculate the transfer matrix for a rectangular barrier

        Parameters
        ----------
        v_out : float
            Wavevector in free space
        v_in : float
            Wavevector in the potential barrier
        d : float
            width of potential barrier

        Returns
        -------
        M : ndarray
        """
        if v_in == 0 or v_out == 0:
            return np.zeros((2, 2))
        m11 = 1j / 2 * (v_out / v_in - v_in / v_out) * sinh(v_in * d) + cosh(v_in * d)
        m12 = -1j / 2 * (v_out / v_in + v_in / v_out) * sinh(v_in * d)
        m21 = 1j / 2 * (v_out / v_in + v_in / v_out) * sinh(v_in * d)
        m22 = 1j / 2 * (v_in / v_out - v_out / v_in) * sinh(v_in * d) + cosh(v_in * d)
        return np.array([[m11, m12], [m21, m22]])

    @staticmethod
    def _transfer_matrix_free(v_out, l):
        """Calculate the transfer matrix for a free space

        Parameters
        ----------
        v_out : float
            Wavevector in free space
        l : float
            Length of the free space

        Returns
        -------
        M : ndarray
        """
        return np.array([[exp(1j * v_out * l), 0], [0, exp(-1j * v_out * l)]])
