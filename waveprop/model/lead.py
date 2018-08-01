# -*- coding: utf-8 -*-
"""
Created on 26 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


This module contains the Lead-object for representing the periodic leads of the model

A infinite periodic lead consists of repeating cells (See model.unit_cell). Due to the
periodicity of the potetnial the blochvector has to be calculated using the transfer-matrix
of the lead-cell.

Notes
-----
The Transmission through the lead is zero for an energy-value resulting in an imaginary bloch-vector

See Also
--------
    - Cell

Examples
--------
Create a new lead from a unit-cell:

    cell = Cell.rectangle_barrier(v=10, a=1.0, d=0.8)
    lead = Lead(cell)

Calculate the value of the transcendental equation, the resulting bloch-vector and the transmission amplitude
for a given energy e:

    e = 10
    lead.set_energy(e)
    trans_eq = lead.transcendent_equation
    k_bloch = lead.bloch_vector
    tm = lead.transfer_matrix
    t = tm.t

"""

import numpy as np
from .unit_cell import Cell
from ..core import TransferMatrix


class Lead:

    def __init__(self, lead_cell=None):
        """ Initialize Lead-instance

        Parameters
        ----------
        lead_cell: Cell, optional
            unit-cell used in the lead
        """
        self._cell = None
        self._e = 0
        self._transfer_matrix = TransferMatrix()
        if lead_cell:
            self.load_cell(lead_cell)

    def load_cell(self, lead_cell):
        """ Load cell used in Lead

        Parameters
        ----------
        lead_cell: Cell
        """
        self._cell = lead_cell

    def load_units(self, approximation_units):
        """ Load units into the lead-cell and store it

        Parameters
        ----------
        approximation_units: array_like of ApproximationUnit
            units to use in the lead_cell
        """
        cell = Cell(approximation_units)
        self._cell = cell

    def set_energy(self, e):
        """ Sets the energy of the lead-model and calculate all values

        Parameters
        ----------
        e : float
            energy of the considered particle moving through the lead
        """
        self._e = e
        self._transfer_matrix.multiple_rectangles(e, self.units)

    @property
    def cell(self):
        """ Cell: the cell-object used for the lead """
        return self._cell

    @property
    def units(self):
        """ array_like of ApproximationUnits: get the unit-objects of the cell used in the lead """
        return self._cell.units

    @property
    def transfer_matrix(self):
        """ TransferMatrix: Calculate and return the transfer-matrix of the cell used in the lead """
        if np.imag(self.bloch_vector) != 0:
            self._transfer_matrix.null()
        return self._transfer_matrix

    @property
    def transcendent_equation(self):
        """ float: get the value of the transcendent equation for the currently set energy """
        return 0.5 * self._transfer_matrix.trace()

    @property
    def bloch_vector(self):
        """ float: get the value of the bloch-vector for the currently set energy"""
        return np.arccos(self.transcendent_equation) / self._cell.a

    def transmission_curve(self, elim, steps=1000):
        """ Calculate the transmission data (e- and t-values) for the cell used in the lead

        Parameters
        ----------
        elim: tuple of float
            energy-range to calculate transmission value
        steps: int, optional
            number of energy-steps, default: 1000

        Returns
        -------
        transmission_data: tuple of array_like
        """
        e_values = np.linspace(*elim, steps)
        t_values = list()
        for e in e_values:
            self.set_energy(e)
            t_values.append(self.transfer_matrix.t)  # calculate transmission  and append to data
        return e_values, t_values

    def potential(self, x0, n, direction=1):
        """ Build the potential data of the lead

        Parameters
        ----------
        x0: float
            start-point of lead
        n: int
            number of cells to build for lead-potential
        direction: int, optional
            positive or negative x-direction of lead, default: positive

        Returns
        -------
        potential-data: tuple of array_like
        """
        a = self._cell.a
        x, y = list(), list()

        # set direction based offset of cell start-point
        x0 = x0-a if direction < 0 else x0

        for _ in range(n):
            # get potential-data of the used cell
            x_cell, y_cell = self._cell.data
            x_cell = np.asarray(x_cell)
            y_cell = np.asarray(y_cell)
            if direction > 0:
                # append next cell on the right
                x = np.append(x, x_cell+x0)
                y = np.append(y, y_cell)
                x0 += a
            else:
                # append next cell on the left
                x = np.append(x_cell+x0, x)
                y = np.append(y_cell, y)
                x0 -= a
        return x, y

    def get_band_structure(self, elim, steps=10000):
        """ Calculate the values of the band-structure of the lead-instance

        Parameters
        ----------
        elim: tuple of float
            energy-range to calculate band-structure
        steps: int, optional
            number of energy-steps, default: 10000

        Returns
        -------
        band_data: tuple of array_like
            e- and k_re- and k_im-values of the band-structure
        """
        # calculate bloch-vector-values
        e_values = np.linspace(*elim, steps)
        k_values = list()
        for e in e_values:
            self.set_energy(e)
            k_values.append(self.bloch_vector)  # calculate bloch-vector and append to data
        k_values = np.asarray(k_values)  # convert to numpy array for mirroring

        # mirror band-structure
        k_values = np.append(-k_values[::-1], k_values)
        e_values = np.append(e_values[::-1], e_values)

        # split values into real- and imaginary-curves
        k_real = np.where(np.imag(k_values) == 0, k_values, np.nan)
        k_imag = np.where(np.imag(k_values) != 0, np.imag(k_values), np.nan)

        return e_values, k_real, k_imag

    def get_transcendent_equation(self, elim, steps=10000):
        """ Calculate the values of the transcendent equation of the lead-instance

        Parameters
        ----------
        elim: tuple of float
            energy-range to calculate the transcendent equation
        steps: int, optional
            number of energy-steps, default: 10000

        Returns
        -------
        equation_values: tuple of array_like
        """
        e_values = np.linspace(*elim, steps)
        eq_values = list()
        for e in e_values:
            self.set_energy(e)
            eq_values.append(self.transcendent_equation)  # calculate transcendent-eq and append to data
        return e_values, eq_values
