# -*- coding: utf-8 -*-
"""
Created on 27 Jun 2018
@author: Dylan Jones

This module contains the Sample-object for representing the sample region of the model

The Sample-object can be loaded using two different methods:
    - sample.load_units:
        loads multiple approximation-units as representation of the samplepotential
    - sample.load_cells:
        loads the approximation units of one or more cell-objects. This can be used for
        fast prototyping of potetnials

See Also
--------
    - Cell

Examples
--------
Create a new sample from a unit-cell:

    cell = Cell.rectangle_barrier(v=10, a=1.0, d=0.8)
    sample = Sample()
    sample.load_cells([cell])

Calculate the value of the transmission amplitude for a given energy e:

    e = 10
    sample.set_energy(e)
    tm = sample.transfer_matrix
    t = tm.t
"""

import copy
import numpy as np
from ..core.transfer_matrix import TransferMatrix


class Sample:

    def __init__(self, approximation_units=None):
        """ Initialize the Sample-instance

        Parameters
        ----------
        approximation_units: array_like of ApproximationUnit, optional
            if not None, load units into model

        """
        self._units = None
        self._width = None
        self._e = 0
        self._transfer_matrix = TransferMatrix()

        if approximation_units:
            self.load_units(approximation_units)

    def load_units(self, approximation_units):
        """ Load units into the sample and store them

        Parameters
        ----------
        approximation_units: array_like of ApproximationUnit
            units to use in sample region
        """
        self._units = approximation_units
        self._width = sum([unit.width for unit in self._units])

    def load_cells(self, cells):
        """ Load cell-objects into the sample region.

        This saves all units used in the cells in the sample-instance

        Parameters
        ----------
        cells: array_like of Cell
            cells used in sample-region
        """
        units = list()
        x0_cell = 0
        for cell in cells:
            cell_units = list(cell.units)
            for unit in cell_units:
                unit_copy = copy.copy(unit)
                unit_copy.x0 += x0_cell
                units.append(unit_copy)
            x0_cell += cell.a

        self.load_units(units)

    def set_energy(self, e):
        """ Sets the energy of the sample-model and calculate all values

        Parameters
        ----------
        e : float
            energy of the considered particle moving through the lead
        """
        self._e = e
        self._transfer_matrix.multiple_rectangles(e, self.units)

    @property
    def xlim(self):
        """ tuple of float: x-range of the sample-potential"""
        x0 = self.units[0].x0
        return x0, x0 + self._width

    @property
    def units(self):
        """ array_like of ApproximationUnit: units of the approximated sample"""
        return self._units

    @property
    def transfer_matrix(self):
        """ TransferMatrix: The trasfermatrix-instance for the sample"""
        return self._transfer_matrix

    def transmission_curve(self, elim, steps=1000):
        """ Calculate the transmission data (e- and t-values) for the sample

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
            t_values.append(self.transfer_matrix.t)
        return e_values, t_values

    def potential(self):
        """ Build the potential data of the sample region

        Returns
        -------
        potential_data: tuple of array_like
        """
        x_data, y_data = list(), list()
        for unit in self._units:
            x, y = unit.build(100)
            x_data = x_data + list(x)
            y_data = y_data + list(y)
        return x_data, y_data
