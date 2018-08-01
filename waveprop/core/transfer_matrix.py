# -*- coding: utf-8 -*-
"""
Created on 17 Jul 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


Contains the TransferMatrix-object for calculating transfer-matrices for different potentials.

The transfer-matrix is implemented for following potential-prototypes:

    - rectangle barrier


Examples
--------

Calculate the transmission for a potential containing two rectangle barriers:

    e = 10  # E_h

    tm1 = TransferMatrix()
    tm1.rectangle(e=e, v=10, d=1.0)

    tm2 = TransferMatrix()
    tm2.rectangle(1=e, v=12, d=1.2)

    total_tm = tm1 * tm2
    t = tm.t

"""

import numpy as np
from numpy import exp, linalg
from cmath import sin, cos
import cmath
from .constants import constants


class TransferMatrix:
    """ transfermatrix object

    object for caclulating the transfermatrix and the transmission-koefficient
    of a scattering region
    """

    def __init__(self, M=None):
        """Create object and initialize transfermatrix"""
        self._M = np.zeros((2, 2)) if M is None else M

    @classmethod
    def null_matrix(cls):
        """ Create zeor-matrix """
        return cls(np.array([[0, 0], [0, 0]]))

    @classmethod
    def one_matrix(cls):
        """ Create unit-matrix """
        return cls(np.array([[1, 0], [0, 1]]))

    def null(self):
        """ Set to an empty transfermatrix """
        self._M = np.zeros((2, 2))

    def ones(self):
        """ Set to unit matrix """
        self._M = np.array([[1, 0], [0, 1]])

    def nan(self):
        """ Set all elements to np.nan """
        self._M = np.array([[np.nan, np.nan], [np.nan, np.nan]])

    @property
    def M(self):
        """ ndarray: transfermatrix of scattering region for the set energy """
        return self._M

    @property
    def t(self):
        """ float: transmission-koefficinet of scattering region for set energy """
        m22 = self._M[1, 1]
        if self.contains_nan:
            print("nan")
            return np.nan
        if m22 == 0:
            return 0
        else:
            return 1 / abs(m22) ** 2

    @property
    def contains_nan(self):
        """ bool: check if transfermatrix contains any np.nan"""
        return np.isnan(self.M).any()

    def trace(self):
        """float: Trace of the transfermatrix"""
        return np.trace(self._M)

    """ ============================== Calculation ============================== """

    def rectangle(self, e, v, d):
        """ Calculate the transfermatrix of a rectangular potential

        Parameters
        ----------
        e : float
            Energy of the particle
        v : float
            Potential strength of the barrier
        d : float
            Width of the barrier
        """
        self._M = self._get_rectangle_tm(e, v, d)

    def multiple_rectangles(self, e, rectangles):
        """ Calculate the transfermatrix for multiple rectangle units

        Parameters
        ----------
        e : float
            energy of the particle
        rectangles : array_like of Rectangle
            list of rectangle-units
        """
        m_total = 1
        for rect in rectangles:
            m = self._get_rectangle_tm(e, rect.v, rect.width)
            m_total = np.dot(m, m_total)
        self._M = m_total

    def diagonalize(self, base_matrix):
        """Represent the current transfermatrix in the base of the given matrix

        Parameters
        ----------
        base_matrix : TransferMatrix
            Base matrix
        """
        if self.contains_nan or base_matrix.contains_nan:
            self.nan()
            return
        q = linalg.eig(base_matrix.M)[1]
        m_diag = np.dot(linalg.inv(q), np.dot(self._M, q))
        self._M = m_diag

    """ ============================================================ """

    def _get_rectangle_tm(self, e, v, d):
        """ Calculate the transfermatrix of a rectangular potential

        Parameters
        ----------
        e : float
            Energy of the particle
        v : float
            Potential strength of the barrier
        d : float
            Width of the barrier

        Returns
        -------
        m: ndarray
            transfer-matrix
        """
        v_out = cmath.sqrt(2 * constants.m * e) / constants.hbar
        # v_out = 1j * v_out if e < 0 else v_out
        if e == v:
            m = np.array([[np.nan, np.nan], [np.nan, np.nan]])
        elif v == 0:
            m = self._transfer_matrix_free(v_out, d)
        else:
            v_in = cmath.sqrt(2 * constants.m * (e - v)) / constants.hbar
            m = self._transfer_matrix_rectangle(v_out, v_in, d)
        return m

    @staticmethod
    def _transfer_matrix_rectangle(v_out, v_in, d):
        """ Calculate the transfer matrix for a rectangular barrier

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
        m11 = 1j / 2 * (v_out/v_in + v_in/v_out) * sin(v_in * d) + cos(v_in * d)
        m12 = 1j / 2 * (v_in/v_out - v_out/v_in) * sin(v_in * d)

        m21 = -1j / 2 * (v_in/v_out - v_out/v_in) * sin(v_in * d)
        m22 = -1j / 2 * (v_out/v_in + v_in/v_out) * sin(v_in * d) + cos(v_in * d)
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

    def __mul__(self, other):
        """ Multiplication dunder-method

        Parameters
        ----------
        other : TransferMatrix

        Returns
        -------
        result: TransferMatrix
        """
        if isinstance(other, TransferMatrix):
            m = np.dot(self.M, other.M)
        else:
            m = other*self.M
        return TransferMatrix(m)

    def __rmul__(self, other):
        """ Reversed multiplication dunder-method"""
        return other.__mul__(self)


