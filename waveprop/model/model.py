# -*- coding: utf-8 -*-
"""
Created on 27 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0

This module contains the main model and methods to load it.

The model can be loaded with either lead-objects, a sample region or both.
If leads are loaded with a sample region, the transfer matrix of the sample will be
diagonalized with the transfer matrix of the leads.

To initialize a model from raw potential data, the method "build_model_from_data" can be used.
It analyzes the data and approximates the components.
Both components can be initialized using cell-objects, too


See Also
--------
    - Cell
    - Lead
    - Sample
    - Approximation

Examples
--------
Initialize model with consisting of only a sample region with two rectangle barriers:

    sample_cells = [sample_cell_1, sample_cell_2]
    model = Model()
    model.load_cell(sample_cells=sample_cells, lead_cell=None)

load the same sample region with leads:

    model.load_cell(sample_cells=sample_cells, lead_cell=lead_cell)

load_approximation units (this is used by the "build_model_from_data"-method:

    model = Model()
    model.load_approximation(sample_approx, cell_approx)

    or

    model = build_model_from_data(data_index=0, spin="up", approx_width=0.25)

"""

import numpy as np
import matplotlib.pyplot as plt
from data import get_data
from .lead import Lead
from .sample import Sample
from ..core import Approximation, TransferMatrix


def build_model_from_data(data_index=0, spin="up", approx_width=0.25):
    """ Builds the model from dataset specified in "data"-folder and -object

    Notes
    -----
    If the Approximation-width is lower then the data-resolution, the width-value will
    be raised to the lowest-resolution

    Parameters
    ----------
    data_index: int, optional
        index of dataset
    spin: str, optional
        spin-value of dataset
    approx_width: float, optional
        width of the approsimation-units

    Returns
    -------
    model: Model
    """
    # load and analyze data (see get_data-method)
    sample_data, lead_cell_data = get_data(data_index, spin)

    # check data-resolution
    sample_step = sample_data[0][1] - sample_data[0][0]
    cell_step = lead_cell_data[0][1] - lead_cell_data[0][0]

    # if the approsimation-width is set lower then the resolution
    # raise the approx.-width to the resolution
    max_approx_width = max(sample_step, cell_step, approx_width)
    if max_approx_width != approx_width:
        print("Changed approximation-width to {:.4f} (sample rate)".format(max_approx_width))
    else:
        print("Using approximation-width of {:.4f}".format(max_approx_width))

    # Approximate the lead- and sample-data
    sample_approx = Approximation(*sample_data)
    sample_approx.rectangles(max_approx_width)
    cell_approx = Approximation(*lead_cell_data)
    cell_approx.rectangles(max_approx_width)

    # initialize the model
    model = Model()
    model.load_approximation(sample_approx, cell_approx)
    return model


class Model:

    def __init__(self):
        """ Initialize the model-instance and declare variables for the lead- and sample-component"""
        self._sample = None
        self._lead = None

    def load_approximation(self, sample_approx=None, lead_cell_approx=None):
        """ Load Approsimation-objects as lead and/or sample

        Parameters
        ----------
        sample_approx: Approximation, optional
            Approsimation-instance of the sample-data
        lead_cell_approx: Approximation, optional
            Approximation-instance of the lead-cell-data
        """
        # if sample_approx is given, load sample
        if sample_approx:
            sample = Sample()
            sample.load_units(sample_approx.units)
            self._sample = sample

        # if lead_cell_approx is given, load lead
        if lead_cell_approx:
            lead = Lead()
            lead.load_units(lead_cell_approx.units)
            self._lead = lead

    def load_cell(self, sample_cells=None, lead_cell=None):
        """ Load cell-objects as lead and/or sample

        Parameters
        ----------
        sample_cells: array_like of Cell or Cell, optional
            list of cells of the sample region.
            will be converted to list if single item
        lead_cell: Cell, optional
            cell of the lead
        """
        # if sample_cells is/are given, load sample
        if sample_cells:
            # convert input to list, if single cell-instance
            sample_cells = [sample_cells] if not isinstance(sample_cells, list) else sample_cells
            sample = Sample()
            sample.load_cells(sample_cells)
            self._sample = sample

        # if lead_cell is given, load lead
        if lead_cell:
            lead = Lead()
            lead.load_cell(lead_cell)
            self._lead = lead

    # def load_lead(self, cell):
    #     lead = Lead()
    #     lead.load_cell(lead_cell)
    #     self._lead = lead
    #
    # def load_sample_cells(self, cells):
    #     sample = Sample()
    #     sample.load_cells(cells)
    #     self._sample = sample

    def set_energy(self, e):
        """ Set the of the model and all sub-components, if loaded

        Parameters
        ----------
        e : float
            energy of the considered particle moving through the model
        """
        # if lead-instance is loaded, set its energy
        if self._lead:
            self._lead.set_energy(e)

        # if sample-instance is loaded, set its energy
        if self._sample:
            self._sample.set_energy(e)

    @property
    def lead(self):
        """ Lead: lead-component of the model"""
        return self._lead

    @property
    def sample(self):
        """ Sample: sample-component of the model"""
        return self._sample

    @property
    def transfer_matrix(self):
        """ Get Transfermatrix of full model

        If lead-component is loaded, the transfer matrix must bediagonalized in the base of the
        lead-transfer -matrix

        Returns
        -------
        tm: TransferMatrix
        """
        lead_tm = self._lead.transfer_matrix if self._lead else TransferMatrix.one_matrix()
        sample_tm = self._sample.transfer_matrix if self._sample else TransferMatrix.one_matrix()

        # diagonalize transfer matrix
        if self._lead:
            # handle bad values
            if lead_tm.t == 0:
                return TransferMatrix.null_matrix()

            # if no sample is loaded, set sample-transfer-matrix to the unitary matrix
            tm_sample = self._sample.transfer_matrix if self._sample else TransferMatrix.one_matrix()
            tm_sample.diagonalize(lead_tm)
        return sample_tm

    def transmission_curve(self, elim, steps=1000):
        """ Calculate the transmission data (e- and t-values) for the model

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

    def potential(self, n_lead_cells=2):
        """ Build the potential data of all components loaded in the model

        Parameters
        ----------
        n_lead_cells: int, optional
            number of cells to plot in each lead, default: 2

        Returns
        -------
        potential-data: tuple of array_like
        """
        sample = self._sample.potential()
        x0, x1 = self._sample.xlim
        if self.lead:
            lead_in = self._lead.potential(x0, n_lead_cells, -1)
            lead_out = self._lead.potential(x1, n_lead_cells, 1)
        else:
            sample_w = abs(x1-x0)
            offset = 0.5 * sample_w
            lead_in = [x0-offset, x0], [0, 0]
            lead_out = [x1, x1+offset], [0, 0]
        x, y = list(), list()
        for _data in [lead_in, sample, lead_out]:
            x = x + list(_data[0])
            y = y + list(_data[1])
        return x, y

    @property
    def potential_lim(self):
        """ tuple of float: Returns the min and max value of the potential of the model"""
        v = self.potential(1)[1]
        return min(v), max(v)

    def get_cell_edges(self, n_lead_cells):
        """ calculates the edges of the cells of each lead

        Parameters
        ----------
        n_lead_cells: int
            number of cells to plot in each lead

        Returns
        -------
        x_edges: array_like
            x-values of the cell-edges
        """
        x0, x1 = self._sample.xlim
        a = self._lead.cell.a
        x = list()
        for i in range(n_lead_cells+1):
            x.append(x0 - i*a)
            x.append(x1 + i*a)
        x.sort()
        return x

    def plot(self, n_lead_cells=5, col="black"):
        """ Plots the potential of the model

        Parameters
        ----------
        n_lead_cells : int, optional
            number of cells to plot in each lead, default: 5
        col: matplotlib.color, optional
            color for plotting, default: black

        Returns
        -------
        fig: matplotlib.figure
            figure-instance of plot
        ax: matplotlib.axis
            axis-object of plot
        """
        fig, ax = plt.subplots()
        x, y = self.potential(n_lead_cells)
        ax.plot(x, y, color=col, lw=1.5)
        if self.lead:
            for x in self.get_cell_edges(n_lead_cells):
                ax.axvline(x, lw=0.5, color="0.5")
        return fig, ax
