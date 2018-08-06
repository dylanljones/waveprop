# -*- coding: utf-8 -*-
"""
Created on 22 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


This module contains a data-container for the raw data of the potential.
At runtime the RawData-instances are created and stored as list.
Analyzed potential-parts can be accessed with the 'get_data'-method

Example
-------
Get spin-up-potential of first dataset:

    index = 0
    spin = 'up'
    sample_data, lead_cell_data = get_data(index, spin=spin)

Attributes
----------
FILE_NAME: str
    .\ + <filename>: must be in same directory as this module
potentials: array_like
    list conatining the RawData-instances from the file specified
"""

import os
import numpy as np
from .potential_analysis import PotentialAnalyzer

#  name of the file containing the potential data
FILE_NAME = r"0.test-VT_AV.dat"  # must be in same directory


class RawData:
    """ Class containing the raw data of the spin-up and -down potential """

    def __init__(self, x, up, down):
        """ Initialization of the data-container

        Parameters
        ----------
        x : array_like
            x-values of the potential
        up : array_like
            y-values for the spin-up potential
        down : array_like
            y-values for the spin-down potential
        """
        self._x = x
        self._up = up
        self._down = down

    @property
    def x(self):
        """ array_like: x-values of the spin-up and -down potentials """
        return self._x

    @property
    def up(self):
        """ y-values for the spin-up potential """
        return self._up

    @property
    def down(self):
        """ y-values for the spin-down potential """
        return self._down

    @property
    def data_up(self):
        """ tuple of array_like: x- and y-values for the spin-up potential """
        return self._x, self._up

    @property
    def data_down(self):
        """ tuple of array_like: x- and y-values for the spin-down potential """
        return self._x, self._down

    @property
    def data(self):
        """ tuple of tuple of array_like: x- and y-values for the spin-up and -down potential """
        return self.data_up, self.data_down

    def x_range(self):
        """ float: width of the potential data """
        return abs(max(self.x)-min(self.x))


def _get_raw_data():
    """ Load the data-file and store the different potentials as RawData-objects

    Returns
    -------
    pots: array_like
        list of RawData-objects containing
    """
    file = os.path.join(os.path.dirname(__file__), FILE_NAME)
    data = np.loadtxt(file).T

    x = data[0]  # first collumn in file is x-data
    v_arrays = data[1:]  # rest are y-values of the potentials
    pots = list()

    pots.append(RawData(x, v_arrays[0], v_arrays[7]))
    pots.append(RawData(x, v_arrays[1], v_arrays[2]))
    pots.append(RawData(x, v_arrays[4], v_arrays[5]))
    # potentials.append(Potential(x, v_arrays[3], v_arrays[6]))
    return pots


def get_data(index, spin="up"):
    """ Get dataset of potential for up- and down-spin and analyze data

    Parameters
    ----------
    index : int
        index of potential
    spin : str, optional
        spin-direction

    Returns
    -------
    sample_data: tuple of array_like
        x- and y-data of the sample-region
    lead_cell_data: tuple of array_like
        x- and y-data of the lead-cell
    """
    global potentials  # potentials = _get_raw_data()

    # get raw-data
    pot = potentials[index]
    if spin.lower() == "up":
        raw_data = pot.data_up
    elif spin.lower() == "down":
        raw_data = pot.data_down
    else:
        return None

    # analyze data and extract sample and lead-cell
    analyzer = PotentialAnalyzer()
    analyzer.load(raw_data)
    sample_data = analyzer.sample_data
    lead_cell_data = analyzer.lead_cell_data

    return sample_data, lead_cell_data


# Read data-file and store potential-data for further use
potentials = _get_raw_data()
