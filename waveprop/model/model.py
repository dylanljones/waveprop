# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones


This module contains the main model of the package, consisting of leads and
an optional sample.

"""
from .kronig_penney import KronigPenney as Lead
from .sample import Sample, OrderedSample, DisorderedSample
from .unitcell import Cell
from waveprop.plotting.plot_utils import plot_model, add_model_plot
from numpy import linalg
import numpy as np


class Model:
    """ Main model of package
    _______________________________________
                |               |
        Lead    |    Sample     |   Lead
    ____________|_______________|__________

    main model consisting of to semi-infinite leads modelled with the
    kronig-penney-model and an optional sample between the two.
    """

    def __init__(self, lead_cell=None, n_bands=3):
        """ Create the model

        Parameters
        ----------
        lead_cell : Cell, optional
            unitcell used in leads
        n_bands : int, optional
            number of bands to calculate
        """
        self._e = 0
        self._lead = None
        self._sample = None

        self.set_lead(lead_cell, n_bands)

    @classmethod
    def from_lead_params(cls, v, a, d):
        """Create model from cell-parameters

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
        model : Model
        """
        cell = Cell(v, a, d)
        return cls(lead_cell=cell)

    """ ================================= Configuration ======================================== """

    def set_lead(self, lead_cell, n_bands):
        """Set leads of the model

        Parameters
        ----------
        lead_cell : Cell
            unitcell used in leads
        n_bands : int
            number of bands to calculate
        """
        self._lead = Lead(lead_cell, n_bands)

    def set_sample(self, cells=None):
        """Set the sample of the model

        Parameters
        ----------
        cells : list of Cell, optional
            unitcells of the scattering region
        """
        if cells is None:
            self._sample = None
        else:
            self._sample = Sample(cells)

    def set_ordered_sample(self, n, cell=None, n_bands=3):
        """Set an ordered Sample region

        Parameters
        ----------
        n : int
            number of unitcells in scattering region
        cell : Cell
            unitcell used in sample
        n_bands : int
            number of bands to calculate
        """
        if isinstance(self._sample, OrderedSample):
            self._sample.new_sample(n, cell)
            self._sample.set_energy(self.e)
        elif cell is not None:
            sample = OrderedSample(cell, n, n_bands)
            self._sample = sample
            self._sample.set_energy(self.e)
        else:
            self._sample = None

    def set_disordered_sample(self, w, n, cell=None):
        """Set a disordered Sample region

        Parameters
        ----------
        w : float
            Disorder strength of sample
        n : int
            number of unitcells in scattering region
        cell : Cell
            unitcell used for calculating random cells in sample region
        """
        if isinstance(self._sample, DisorderedSample):
            self._sample.new_sample(w, n, cell)
            self._sample.set_energy(self.e)
        else:
            if cell is None:
                cell = self.lead_cell
            sample = DisorderedSample(cell, w, n)
            self._sample = sample
            self._sample.set_energy(self.e)

    def is_ordered(self):
        """Check if sample is an ordered sample

        Returns
        -------
        is_ordered : bool
        """
        return True if isinstance(self.sample, OrderedSample) else False

    def is_disordered(self):
        """Check if sample is a disordered sample

        Returns
        -------
        is_disordered : bool
        """
        return True if isinstance(self.sample, DisorderedSample) else False

    def set_energy(self, e):
        """Set energy of particle moving through Model

        Sets the energy of leads (kronig-penney-model) and the sample region

        Parameters
        ----------
        e : float
            Energy of particle
        """
        self.e = e
        self._lead.set_energy(self.e)
        if self._sample is not None:
            self._sample.set_energy(self.e)

    """ ==================================== Properties ======================================== """

    @property
    def e(self):
        """float: Currently set energy"""
        return self._e

    @e.setter
    def e(self, e):
        if e == np.nan:
            e = 0
        self._e = abs(e)

    @property
    def lead(self):
        """KronigPenney: Lead object used in model"""
        return self._lead

    @property
    def sample(self):
        """Sample: Sample object used in model"""
        return self._sample

    @property
    def lead_cell(self):
        """Cell: unitcell used in leads"""
        return self._lead.cell

    @property
    def sample_cell(self):
        """Cell: unitcell used in sample, if ordered"""
        if isinstance(self.sample, OrderedSample):
            return self._sample.cell

    @property
    def bands(self):
        """list of list: energy bands of the leads"""
        return self.lead.bands

    @property
    def gaps(self):
        """list of list: energy gaps of the leads"""
        return self._lead.gaps

    @property
    def sample_bands(self):
        """list of list: energy bands of the sample, if ordered"""
        if isinstance(self.sample, OrderedSample):
            return self.sample.bands
        else:
            return None

    @property
    def sample_gaps(self):
        """list of list: energy gaps of the sample, if ordered"""
        if isinstance(self.sample, OrderedSample):
            return self.sample.gaps
        else:
            return None

    @property
    def sample_length(self):
        """float: total length of the sample"""
        if self._sample is None:
            return 0
        else:
            return self._sample.length

    @property
    def k(self):
        """float: blochvector of model"""
        return self._lead.k

    @property
    def sample_k(self):
        """float: blochvector of sample, if ordered"""
        if isinstance(self.sample, OrderedSample):
            k = 0
            if np.imag(self.k) == 0:
                k = self.sample.k
            return k
        else:
            return None

    @property
    def t(self):
        """float: Transmission koefficient of model"""
        return self.transmission()

    @property
    def M(self):
        """ndarray: Transfermatrix of model"""
        return self._transfer_matrix()

    """ =================================== Main Methods ======================================= """

    def transmission(self, e=None):
        """Calculate the transmission koefficient for a given energy e

        Parameters
        ----------
        e : float
            energy of particle

        Returns
        -------
        t : float
            transmission koefficient
        """
        if e is None:
            e = self.e
        self.set_energy(e)
        if np.imag(self.k) != 0:
            return 0
        m = self._transfer_matrix()
        m22 = m[1, 1]
        if m22 == 0:
            return np.nan
        else:
            return 1 / (np.abs(m22) ** 2)

    def transmission_curve(self, xlim=None, n=2000):
        """Calculate transmission values for energies in a given range

        Parameters
        ----------
        xlim :list of float
            energy range for calculating transmission curve, consisting of
            the start and end value.
        n : int
            number of energy levels to calculate

        Returns
        -------
        data : list of ndarray
            e and t data of the transmission curve
        """
        if xlim is None:
            xlim = [0, 30]
        e_values = np.linspace(*xlim, n)
        t_values = []
        for e in e_values:
            t_values.append(self.transmission(e))
        return e_values, t_values

    def get_band(self, i, offset=None):
        """Get the start and end of the i-th energy band

        Parameters
        ----------
        i : int
            index of band
        offset : float
            scaling factor for different representations

        Returns
        -------
        band_energies :obj:"list" of float
        """
        e0, e1 = self.bands[i - 1]
        if offset:
            delta = self.get_band_size(i) * offset
            e_range = [e0 - delta, e1 + delta]
        else:
            e_range = [e0, e1]
        return e_range

    def get_band_size(self, i):
        """Get the size of the i-th energy band

        Parameters
        ----------
        i : int
            index of band

        Returns
        -------
        size : float
        """
        band = self.get_band(i)
        return abs(band[1]-band[0])

    def in_sample_band(self, e):
        """Check if given energy e is in energy band of sample

        Parameters
        ----------
        e : float
            energy level

        Returns
        -------
        in_band : bool
        """
        if self.sample.type() != OrderedSample:
            return False
        for band in self.sample_bands:
            in_band = (band[0] <= e) and (e <= band[1])
            if in_band:
                return True
        return False

    def in_lead_band(self, e):
        """Check if given energy e is in energy band of leads

        Parameters
        ----------
        e : float
            energy level

        Returns
        -------
        in_band : bool
        """
        for band in self.bands:
            if self.in_band(e, band):
                return True
        return False

    @staticmethod
    def in_band(e, band):
        """Check if given energy e is in given band

        Parameters
        ----------
        e : float
            energy level
        band : array_like of float
            band to check

        Returns
        -------
        in_band : bool
        """
        if (band[0] <= e) and (e <= band[1]):
            return True
        else:
            return False

    def total_bands(self):
        """Calculate total band structure

        If sample is an ordered sample, the band structure is an overlap of the bands
        of the leads and the bands of the sample region

        Returns
        -------
        total_bands : list of list
        """
        if not self.is_ordered():
            return self.bands
        total_bands = []
        for band in self.bands:
            for sample_band in self.sample_bands:
                e0, e1 = sample_band
                lower_in = self.in_band(e0, band)
                upper_in = self.in_band(e1, band)
                both_in = lower_in and upper_in
                if both_in:
                    new_band = [e0, e1]
                    total_bands.append(new_band)
                elif lower_in and not upper_in:
                    new_band = [e0, band[1]]
                    total_bands.append(new_band)
                elif not lower_in and upper_in:
                    new_band = [band[0], e1]
                    total_bands.append(new_band)
        return total_bands

    """ ====================================== Other =========================================== """

    def _transfer_matrix(self):
        """Calculate the total transfermatrix of the system

        First the eigenbasis q of the transfermatrix of the leads is calculated.
        To get the total matrix, the transfermatrix of the sample is represented in
        the base of the lead-transfermatrix

        Returns
        -------
        M : ndarray
            transfermatrix of the whole system
        """
        q = linalg.eig(self.lead.M)[1]
        if self._sample is None:
            m_s = 1
        else:
            m_s = self._sample.M
        return np.dot(linalg.inv(q), np.dot(m_s, q))

    def build(self, cell_resolution, n_cells=None):
        """Build potential data for plotting the model

        Parameters
        ----------
        cell_resolution : int
            Number of points in each cell
        n_cells : int
            Number of lead-cells in build

        Returns
        -------
        potential_data : list of list
            potential data of the two leads and the sample
        """
        n_sample = 2
        if self.sample is not None:
            n_sample = self.sample.n
        if n_cells is None:
            n_cells = max(2, int(n_sample / 25))
        lead = self.lead
        a = lead.cell.a

        x0 = - n_cells * a
        lead1 = self.lead.build(x0, n_cells, cell_resolution)

        if self.sample is not None:
            sample = self.sample.build(cell_resolution=cell_resolution)
        else:
            sample = [np.asarray([]), np.asarray([])]
        x0 = self.sample_length
        lead2 = self.lead.build(x0, n_cells, cell_resolution)
        return lead1, lead2, sample

    def show(self, font=None, cell_resolution=None, n_cells=None):
        """ Creates plot of current model

        Parameters
        ----------
        font: dict
            font for plot
        cell_resolution : int
            Number of points in each cell
        n_cells : int
            Number of lead-cells in build

        Returns
        -------
        plot : matplotlib figure
        """
        plt = plot_model(self, font, cell_resolution, n_cells)
        return plt

    def add_plot(self, plt, cell_res, n_cells=None, label=None, sample_color=None):
        _, ylim = add_model_plot(plt.ax, self, cell_res, n_cells, label, sample_color)
        return ylim

    def __str__(self):
        out = "\nLEAD:\n"
        out += str(self.lead)
        out += "SAMPLE: "
        out += str(self.sample) + "\n"
        out += "________________________"
        return out
