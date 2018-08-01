# -*- coding: utf-8 -*-
"""
Created on 27 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


contains the Approximation-object for calculating and storing approximation of the potential data

Available approximation-types:

    - Rectangle approximation

Examples
--------
Approximate sinus-function with rectangles and build approximated data:

    x = np.linspace(0, 2*np.pi, 1000)
    y = np.sin(x)

    approx = Approximation(x, y)
    unit_width = np.pi/10
    approx.rectangles(width=unit_width, sample_pos="c")
    x_approx, y_approx = approx.build()
"""

import matplotlib.pyplot as plt
from .units import Rectangle


class Approximation:

    def __init__(self, x, y):
        """  Initialize the Approximation-instance

        Parameters
        ----------
        x : array_like
            x-values of the data to approximate
        y : array_like
            y-values of the data to approximate
        """
        self._x = x
        self._y = y
        self._units = list()  #: list of units in approximation

    @property
    def units(self):
        """ array_like: unit-objects used in approximation """
        return self._units

    @property
    def unit_width(self):
        """ float: width of the units"""
        return self.units[0].width

    def rectangles(self, width=None, y0=0, sample_pos="c"):
        """ Approximate data using the rectangle-unit

        Parameters
        ----------
        width : float
            width of the units for the approximation
        y0 : float, optional
            y-offset, default: 0
        sample_pos: str, optional
            position of sampling value in unit:
            l: left, c: center, r: right
            default: center

        Returns
        -------
        units: array_like
        """
        self._units = list()

        # calculate step-size for units
        n_data_points = len(self._x)
        i_steps = self._get_index_steps(width, 0.01)
        step_size = self._x[i_steps]-self._x[0]

        # get offset corresponding to sampling-position
        offset = 0
        if sample_pos == "c":
            offset = i_steps/2
        elif sample_pos == "r":
            offset = i_steps - 1
        elif sample_pos == "l":
            offset = 0

        # calculate approximation
        if i_steps > 0:
            i = 0
            while i < n_data_points - offset:
                i_pos = i
                x0 = self._x[i_pos]  # unit position
                i_s = int(i_pos + offset)  # index of sampling value
                sample_val = self._y[i_s]  # sample value for unit
                rect = Rectangle(x0, step_size, sample_val-y0)  # create Rectangle-unit
                self._units.append(rect)
                i += i_steps

        return self._units

    def build(self, steps_per_unit=100, draw_to_floor=False):
        """ Build data of the calculated approximation

        Parameters
        ----------
        steps_per_unit : int, optional
            data points per unit, default: 100
        draw_to_floor : bool, optional
            if true, draw edges of units to zero
            default: False

        Returns
        -------
        data: tuple of array_like
        """
        x_values = list()
        y_values = list()
        for unit in self._units:
            x, y = unit.build(steps_per_unit, draw_to_floor)
            x_values += list(x)
            y_values += list(y)
        return x_values, y_values

    def show(self):
        """ Plot the original data and the approximation"""
        fig, ax = plt.subplots()
        ax.plot(self._x, self._y, color="r", lw=1, label="Original data")
        ax.plot(*self.build(), color="black", lw=0.5, label="Approximation")
        ax.legend()
        plt.show()

    def _get_index_steps(self, width, default_percent):
        """ Calculate the number of indices per unit

        Parameters
        ----------
        width : float
            width of the units
        default_percent : float
            if width=None, use percentage of number of data-points
        Returns
        -------
        steps: int
        """
        data_step_size = abs(self._x[2]-self._x[1])
        width = default_percent * len(self._x) if width is None else width
        i_steps = int(width/data_step_size)
        i_steps = 1 if i_steps < 1 else i_steps
        return i_steps
