# -*- coding: utf-8 -*-
"""
Created on 24 Mar 2018
@author: Dylan Jones

This module holds the Curve-object and helper methods

"""
import numpy as np
from scipy.optimize import curve_fit

pm_str = chr(177)  # unicode symbol for plus-minus


def function_fit(func, x, y, p0=None, maxfev=800):
    """Fits data to the given function and calculates errors

    Parameters
    ----------
    func : method
        function to fit data to
    x : ndarray
        x-data for fitting
    y : ndarray
        y-data for fitting
    p0 : array_like of float, optional
        starting values for fit parameters
    maxfev : int, optional
        maximum iterations

    Returns
    -------
    popt : ndarray
        list of optimized fit parameters
    errs : ndarray
        list of error-values for fit parameters55
    """
    popt, pcov = curve_fit(func, x, y, p0=p0, maxfev=maxfev)
    errs = np.sqrt(np.diag(pcov))
    return popt, errs


def split_curves(curves, param):
    """Split curves into lists sorted by the given parameter

    Parameters
    ----------
    curves : array_like of Curve
        all curves
    param : str
        key of parameter to sort after

    Returns
    -------
    split_curves : list of list of Curve
    """
    out = {}
    for curve in curves:
        val = curve.params[param]
        if val not in out.keys():
            out.update({val: [curve]})
        else:
            out[val].append(curve)
    return list(out.values())


def filter_curves(curves, param, values, exlude_vals=None):
    """Filter curves after given parameter and/or exclude values

    Parameters
    ----------
    curves : array_like of curves
        all curves
    param : str
        key of parameter to use for filtering
    values : array_like of float
        list fo values to include
    exlude_vals : array_like of float
        list fo values to exclude
    Returns
    -------
    filtered_curves : list of Curve
    """
    if not param:
        return None
    out = []
    for curve in curves:
        val = curve.params.get(param)
        in_range = (values[0] <= val) and (val <= values[1])
        if in_range:
            if exlude_vals:
                if val in exlude_vals:
                    continue
            out.append(curve)

    return out


def get_curves(curves, param, values):
    """Get curves with values for given parameter

    Parameters
    ----------
    curves : array_like of curves
        all curves
    param : str
        key of parameter to use for filtering
    values : array_like of float
        list fo values to include
    Returns
    -------
    filtered_curves : list of Curve
    """
    out = []
    for curve in curves:
        val = curve.params.get(param)
        if val in values:
            out.append(curve)
    return out


class Curve:
    """ Curve object for holding plotting-data.

    Curve object with data points, errors and fit data.
    """

    def __init__(self, params=None, func=None):
        """Create curve instance

        Parameters
        ----------
        params : dict, optional
            Dictionary of curve parameters
        func : method, optional
            method for fitting data
        """
        self._x = []
        self._y = []
        self._y_errs = []
        self._x_errs = []

        self._func = func
        self._fit = None
        self._fit_Values = [[], []]

        self._params = params
        self._label = ""

    def add_point(self, x, y, y_err=None, x_err=None):
        """Add single data-point

        Parameters
        ----------
        x : float
        y : float
        y_err : float, optional
        x_err : float, optional
        """
        self._x.append(x)
        if x_err is not None:
            self._x_errs.append(x_err)
        self._y.append(y)
        if y_err is not None:
            self._y_errs.append(y_err)

    def add_values(self, x_values, y_values, y_errs=None, x_errs=None):
        """Add multiple data points

        Parameters
        ----------
        x_values : array_like of float
        y_values : array_like of float
        y_errs : array_like of float, optional
        x_errs : array_like of float, optional
        """
        self._x += list(x_values)
        if x_errs:
            self._x_errs += list(x_errs)
        self._y += list(y_values)
        if y_errs:
            self._y_errs += list(y_errs)

    def f(self, x):
        """Get value of fit-function for given x

        Parameters
        ----------
        x : float

        Returns
        -------
        y : float
        """
        return self.func(x, *self.popt)

    """ ======================================================================================== """

    @property
    def params(self):
        """dict : Parameters of curve"""
        return self._params

    @params.setter
    def params(self, params):
        if isinstance(params, dict):
            self._params = params

    """ data properties """

    @property
    def x(self):
        """ndarray : x-values of curve"""
        return np.asarray(self._x)

    @property
    def y(self):
        """ndarray : y-values of curve"""
        return np.asarray(self._y)

    @property
    def x_errs(self):
        """ndarray : x-errors of curve"""
        return np.asarray(self._x_errs)

    @property
    def y_errs(self):
        """ndarray : y-errors of curve"""
        return np.asarray(self._y_errs)

    @property
    def data(self):
        """ tuple of tuple of ndarray: data of curve"""
        points = self.x, self.y
        errs = self.x_errs, self.y_errs
        return points, errs

    @property
    def points(self):
        """list of tuple : data points of curve"""
        return zip(self.x, self.y)

    @property
    def limits(self):
        """tuple of tuple of float : limits of curve"""
        xlim = min(self.x), max(self.x)
        ylim = min(self.y), max(self.y)
        return xlim, ylim

    """ fit properties """

    @property
    def popt(self):
        """list of float : Optional parameters of curve fit"""
        return self._fit[0]

    @property
    def errs(self):
        """list fo float : Errors of parameters of fit"""
        return self._fit[1]

    @property
    def fit_data(self):
        """tuple of ndarray : x and y data of built fit"""
        return self.fit_x, self.fit_y

    @property
    def fit_x(self):
        """ndarray : x-values of fit"""
        return self._fit_Values[0]

    @property
    def fit_y(self):
        """ndarray : y-values of fit"""
        return self._fit_Values[1]

    @property
    def fit_error(self):
        """ndarray : difference between data points and fit data"""
        return self._fit_error()

    @property
    def has_fit(self):
        """bool : True, if curve has fit-function"""
        if not self._fit:
            return False
        else:
            return True

    @property
    def n(self):
        """int : Number of data points"""
        return len(self._x)

    """ other properties """

    @property
    def error_limits(self):
        """tuple of float : Limits of fit error"""
        return min(self.fit_error), max(self.fit_error)

    @property
    def label(self):
        """str : label of curve"""
        return self._label

    @property
    def func(self):
        """method : currently set fit-function"""
        return self._func

    """ ======================================================================================== """

    def find_max(self):
        """Find maximum of data

        Returns
        -------
        max_point : tuple of float
            Data of maximum
        """
        y_max = max(self._y)
        i_max = self._y.index(y_max)
        return self._x[i_max], y_max

    def find_min(self):
        """Find minimum of data

        Returns
        -------
        min_point : tuple of float
            Data of minimum
        """
        y_max = min(self._y)
        i_max = self._y.index(y_max)
        return self._x[i_max], y_max

    def get_values(self, xlim):
        """Get curve data in given x-range

        Parameters
        ----------
        xlim : array_like of float
            x-range of data
        Returns
        -------
        data : tuple of ndarray
            data in x-range
        """
        out_x, out_y = [], []
        for x, y in zip(self._x, self._y):
            in_range = (xlim[0] <= x) and (x <= xlim[1])
            if in_range:
                out_x.append(x)
                out_y.append(y)
        return np.asarray(out_x), np.asarray(out_y)

    """ ======================================================================================== """

    def set_fit_function(self, func):
        """Set function for fitting curve-data

        Parameters
        ----------
        func : method
        """
        self._func = func

    def fit(self, func=None, p0=None, offset=None, limits=None, accuracy=3, label_errs=False, x_range=None):
        """fit curve data, set optimal parameters and errors and build fit data

        Parameters
        ----------
        func : method, optional
            function for fitting data, if not given use preiosly set function
        p0 : array_like of float, optional
            initial fit-parameters
        offset : float, optional
            offset of building fit
        limits : array_like of float, optional
            limits of building fit
        accuracy : int
            decimal points of fit parameters in label
        label_errs : bool, optional
            set curve label to fit errors
        x_range : array_like of float, optional
            fit only with data in thsi range
        """
        if func is not None:
            self.set_fit_function(func)

        if x_range is not None:
            x, y = self.get_values(x_range)
        else:
            x, y = self.x, self.y

        try:
            popt, errs = function_fit(self._func, x, y, p0=p0)
        except Exception as e:
            print(e)
        else:
            if label_errs:
                self._build_label(popt, errs, accuracy)
            else:
                self._build_label(popt, acc=accuracy)
            self._build_fit(popt, errs, offset=offset, limits=limits)

    def popt_str(self, i, name=None, acc=3, latex=False):
        """Build string of optimal parameters

        Parameters
        ----------
        i : int
            index of fit-parameter
        name : str, optional
            name of parameter
        acc : int
            decimal points of fit parameter
        latex : bool, optional
            build latex string, if true

        Returns
        -------
        popt_str: str
        """
        param = self.popt[i]
        err = self.errs[i]
        string = name + " = " if name else ""
        string += "{:.{acc}f} ".format(param, acc=acc)
        string += "\pm" if latex else pm_str
        string += " {:.{acc}f}".format(err, acc=acc)
        string = "$ " + string + " $" if latex else string
        return string

    """ ======================================================================================== """

    def rebuild_fit(self, x):
        """Rebuild fit-data for given x-values

        Parameters
        ----------
        x : ndarray
        """
        y = self._func(x, *self.popt)
        self._fit = [self.popt, self.errs]
        self._fit_Values = [x, y]

    def _build_fit(self, popt, errs, offset=None, limits=None, n=1000):
        """Helper method for builds fit-data from optimal fit parameters

        Parameters
        ----------
        popt : array_like of float
            optimal fit parameters
        errs : array_like of float
            errors of fit parameters
        offset : float, optional
            offset of building fit
        limits : tuple of float, optional
            limits of building fit
        n : int
            number of points to build
        """
        if offset is None and not limits:
            offset = 0.3

        with_limits = (offset is None) and limits
        with_offset = (offset is not None) and (not limits)

        if with_offset:
            x0, x1 = min(self.x) - offset, max(self.x) + offset
        elif with_limits:
            x0, x1 = limits[0], limits[1]
        else:
            x0, x1 = min(self.x), max(self.x)

        x = np.linspace(x0, x1, n)
        y = self._func(x, *popt)
        self._fit = [popt, errs]
        self._fit_Values = [x, y]

    def _build_label(self, popt, errs=None, acc=3):
        """Helper method to build label from fit parameters

        Parameters
        ----------
          popt : array_like of float
            optimal fit parameters
        errs : array_like of float
            errors of fit parameters
        acc : int
            decimal points of fit parameter
        """
        label = "  Fit: "
        if errs:
            for p, e in zip(popt, errs):
                label += "{:{len}.{acc}f}".format(p, len=acc + 2, acc=acc)
                label += " $\pm$ {:{length}.{acc}f}, ".format(e, length=acc + 3, acc=acc)
        else:
            for p in popt:
                label += "{:{length}.{acc}f}, ".format(p, length=acc + 3, acc=acc)
        self._label += label[:-2]

    def _fit_error(self):
        """Helper method for calculating differenc between data and fit function

        Returns
        -------
        difference : ndarray
        """
        return self._func(self.x, *self.popt) - self.y

    def __str__(self):
        out = "Curve-Object ({} Points)\n".format(len(self.x))
        if self._params:
            out += " -Parameters:\n"
            for p, val in self._params.items():
                out += "   -{}: {:.2f}\n".format(p, val)

        if self.has_fit:
            for i in range(len(self.popt)):
                out += "  -p_{}: ".format(i) + self.popt_str(i) + "\n"
        return out
