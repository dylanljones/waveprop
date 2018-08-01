# -*- coding: utf-8 -*-
"""
Created on 26 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


This module contains the PotentialAnalyzer-object to find the sample- and leadcell. It uses the
LeadFitter-object to assure smooth transissions between the building blocks of the potential.

Examples
--------
Find the sample and lead-data from raw potential values:

    analyzer = PotentialAnalyzer()
    analyzer.load(raw_data)
    sample_data = analyzer.sample_data
    lead_cell_data = analyzer.lead_cell_data
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema


class PotentialAnalyzer:

    """ Analyzer-object to find data of sample-region and lead-cell

    Finds the edges of the sample region by compairing the minimas of the potential.
    The data-region of the leads is fitted with an oszillating function to assure smooth
    transissions between the individual cells and the sample-data with the leads.
    """

    def __init__(self):
        """ Initialize attributes and objects of instance"""
        self._raw_data = None

        self._lead_fit = LeadFitter()  #: class to fit lead-data
        self._i_mins = None  #: indices of potential minimas
        self._sample_edges = None  #: x- and y-values of the edges of the sample

        self._sample_data = None  #: data of the sample
        self._lead_cell_data = None  #: data of the lead-cell

    def load(self, data):
        """ Load data and run analysis

        Parameters
        ----------
        data : tuple of array_like
            x- and y-values of potential
        """
        self._raw_data = data
        self._analyze()

    @property
    def x(self):
        """ array_like: x-values of the potential data"""
        return self._raw_data[0]

    @property
    def y(self):
        """ array_like: y-values of the potential data"""
        return self._raw_data[1]

    @property
    def sample_data(self):
        """ tuple of array_like: x- and y-values of the sample """
        return self._sample_data

    @property
    def lead_cell_data(self):
        """ tuple of array_like: x- and y-values of the lead-cell """
        return self._lead_cell_data

    def _analyze(self):
        """ Analyze data and find sample- and lead-cell-data

        Calculate indices of minimas for further use
        """
        self._i_mins = argrelextrema(self._raw_data[1], np.less)[0]
        self._get_sample()
        self._get_lead()  # must be called after _get_sample!

    def _get_sample(self):
        """ Find data of the sample and store it """

        # find edges of the sample and store data
        i0, i1 = self._get_sample_edges()
        x_sample = self.x[i0:i1+1]
        y_sample = self.y[i0:i1+1]
        self._sample_data = x_sample, y_sample

        # get the start- and end-point of the sample data and store them
        start_point = x_sample[0], y_sample[0]
        end_point = x_sample[-1], y_sample[-1]
        self._sample_edges = start_point, end_point

    def _get_sample_edges(self, amp_error=0.02):
        """ Calculate the indices of the edges of the sample-region

        Compares the minimas of the potential-data.
        At the beginning of the sample, the minima of the potential should vary

        Parameters
        ----------
        amp_error: float, optional

        Returns
        -------
        egde-indices: tuple of int
        """
        amp = abs(max(self.y)-min(self.y))*0.5
        thresh = amp_error * amp

        lead_min = self.y[self._i_mins[0]]
        _prev = None
        i0 = i1 = None
        for i in self._i_mins:
            err = abs(lead_min - self.y[i])
            if i0 is None and err > thresh:
                i0 = _prev
            if i0 is not None and err < thresh:
                i1 = i
                break
            _prev = i
        return i0, i1

    def _get_lead(self):
        """ Find data of the the lead unit-cell and store it

        Note
        ----
        _get_sample must be called first!

        For achieving clean transissions between the unit-cells, the data is fitted with an
        oszillating function (sin) using the LeadFitter instance.
        """

        # Get data-point of the beginning of the sample.
        # Since a fit-function with a center at the minimum value is used, this can be used
        # to shift the fit-function for a clean transmission between the lead and sample
        x0, y0 = self._sample_edges[0]

        # only use first quarter of data for fitting the lead
        data_to_fit = self._data_range(0, 0.25)

        self._lead_fit.set_data(data_to_fit)
        popt = self._lead_fit.popt

        # store the period and x-offset from the optimal parameters
        p, cell_xoffset = popt[1], popt[2]

        # to fix fit-error, pull the data to the sample-startpoint
        y_err = y0 - self._lead_fit.f(x0)
        x_offset = cell_xoffset

        # build fit-result and store lead cell-data
        x, y = self._lead_fit.build([0, p], x_offset=x_offset, y_offset=y_err)
        self._lead_cell_data = x, y

    def _data_range(self, f0, f1):
        """ Get data in given range

        Parameters
        ----------
        f0: float
            start-fraction of data
        f1: float
            end-fraction of data
        Returns
        -------
        data_roi: tuple of array_like
        """

        n_values = len(self._raw_data[0])
        i0 = int(f0 * n_values)
        i1 = int(f1 * n_values)

        x = self._raw_data[0][i0:i1]
        y = self._raw_data[1][i0:i1]
        data_roi = x, y

        # x_range = min(self.x)*f0, max(self.x)*f1
        # data_roi = self._trunc_data(*x_range)

        return data_roi

    # def _trunc_data(self, x0, x1):
    #     """ Truncate data
    #
    #     Parameters
    #     ----------
    #     x0 : float
    #     x1 : float
    #
    #     Returns
    #     -------
    #     data_roi: tuple of array_like
    #     """
    #     i0 = np.abs(self._raw_data - x0).argmin()
    #     i1 = np.abs(self._raw_data - x1).argmin()
    #     x = self._raw_data[0][i0:i1]
    #     y = self._raw_data[1][i0:i1]
    #     return x, y


class LeadFitter:

    def __init__(self):
        """ Initialize the LeadFitter-instance """
        self._data = None
        self._popt = None  #: optimal fit-parameters
        self._errs = None  #: errors of the fit-parameters

    @property
    def popt(self):
        """ array-like: optimal parameters of curve-fit"""
        return self._popt

    def set_data(self, data):
        """ Set data and calculate fit

        Parameters
        ----------
        data : tuple of array_like
            x- and y-values of data
        """
        self._data = data
        self._fit()

    def build(self, xlim, steps=1000, x_offset=0, y_offset=0):
        """ Build x- and y-data of the fit-result

        Parameters
        ----------
        xlim : array_like
            x-range for building fit
        steps : int, optional
            number of x-steps
            default: 1000
        x_offset : float, optional
        y_offset : float, optional

        Returns
        -------
        fit_data: tuple of array_like
            x- and y-values of fit
        """
        x = np.linspace(*xlim, steps)
        y = self.f(x + x_offset) + y_offset
        return x, y

    def f(self, x, x0=None, y0=None):
        """ fit-function with optimal parameters

        Parameters
        ----------
        x : float or array_like
            x-argument of fit-function
        x0 : float, optional
            x-offset
        y0 : float, optional
            y-offset

        Returns
        -------
        y : float
        """
        _a, _p, _x0, _y0 = self._popt
        _x0 = x0 if x0 else _x0
        _y0 = y0 if y0 else _y0
        return self._fit_func(x, _a, _p, _x0, _y0)

    @staticmethod
    def _fit_func(x, a, p, x0, y0):
        """ Fit-function for approximating lead-cell

        Parameters
        ----------
        x : float or array_like
            x-argument
        a : float
            amplitude of oszillation
        p : float
            period of oscillation
        x0 : float
            x-offset of function
        y0 : float
            y-offset of function

        Returns
        -------
        y : float
        """
        return a * (0.5 + np.sin(2*np.pi/p * (x-x0) - np.pi/2)) + y0

    def _fit(self):
        """ Fit data with fit-function and store resulting optional parameters and errors """
        p0 = self._estimate_parameters(self._data)
        self._popt, pcov = curve_fit(self._fit_func, *self._data, p0=p0, maxfev=2000)
        self._errs = np.sqrt(np.diag(pcov))

    @staticmethod
    def _estimate_parameters(data):
        """ Estimate initial values of fit parameters

        Parameters
        ----------
        data : tuple of array_like
            x- and y-values of data

        Returns
        -------
        p0: array_like
        """
        x, y = data
        mean = sum(y) / len(y)
        i_mins = argrelextrema(y, np.less)[0]
        min_x = [x[i] for i in i_mins]
        diffs = [min_x[i + 1] - min_x[i] for i in range(len(min_x) - 1)]

        a0 = abs(max(y) - min(y)) / 2
        p0 = sum(diffs) / len(diffs)
        x0 = x[i_mins[0]]
        y0 = mean-a0/2
        return [a0, p0, x0, y0]
