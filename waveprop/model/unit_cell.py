# -*- coding: utf-8 -*-
"""
Created on 28 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


Contains Cell-object for use in lead or sample

Can be constructed using approxiation-units or values for a rectangle barrier

Examples
--------
Create a cell consisting of a rectangle barrier and some free space:

    cell = Cell.rectangle_barrier(v=10, a=1.0, d=0.8)

Create a cell from approximation-units:

    cell = Cell(units)

"""
from ..core import Rectangle


class Cell:

    def __init__(self, approximation_units=None):
        """ Initialize Cell-object

        Parameters
        ----------
        approximation_units: array_like of ApproximationUnit, optional
            if not None, load unit-objects approximating cell
        """
        self._units = None
        self._a = None
        if approximation_units:
            self.load_units(approximation_units)

    @classmethod
    def rectangle_barrier(cls, v, a, d):
        """ Create a unitcell with a rectangular barrier

        This constructor can be used to quickly build a cell the size of a
        with a rectangular barrier with with d and height v.

        Parameters
        ----------
        v: float
            height of the barrier
        a: float
            width of the cell
        d : float
            width of the barrier in the cell

        Returns
        -------
        unitcell: Cell
        """
        barrier = Rectangle(0, d, v)  #: Create rectangle-unit for the barrier in the cell
        free = Rectangle(d, a-d, 0)  #: Create rectangle-unit for the free space in the cell
        cell = cls()
        cell._units = [barrier, free]
        cell._a = a
        return cell

    def load_units(self, approximation_units):
        """ Load approximation units into the cell-instance

        Parameters
        ----------
        approximation_units: array_like of ApproximationUnit
            unit-objects approximating cell
        """
        self._units = approximation_units
        self._a = sum([unit.width for unit in self._units])

    @property
    def data(self):
        """ tuple of array_like: x- and y-values of the approximated cell"""
        return list(self._get_unit_data())

    @property
    def a(self):
        """ float: width of the cell"""
        return self._a

    @property
    def units(self):
        """ array_like: approximation-units of the cell"""
        return self._units

    def _get_unit_data(self):
        """ Builds the approximation data of the cell

        Returns
        -------
        data: tuple of array_like
            x- and y-values of the approximated cell
        """
        x_data, y_data = list(), list()
        for unit in self._units:
            x, y = unit.build(100)
            x_data = x_data + list(x)
            y_data = y_data + list(y)
        return x_data, y_data
