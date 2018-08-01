# -*- coding: utf-8 -*-
"""
Created on 28 Jun 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0


This module contains the Unit-objects for approximating the potential data and
the unit-prototype class.

Available units:

    - Rectangle

"""

import numpy as np


class ApproximationUnit:
    """ Base Matrix for Approximation-unit used in the potential-approximation"""

    def __init__(self, x0, width):
        """ Initialize the Unit-instance

        Parameters
        ----------
        x0 : float
            position of unit
        width : float
            width of unit
        """
        self._x0 = x0
        self._width = width

    @property
    def x0(self):
        """ float: x-value of the start of the unit """
        return self._x0

    @x0.setter
    def x0(self, value):
        self._x0 = value

    @property
    def width(self):
        """ float: width of the unit"""
        return self._width

    def build(self, steps, zero_edge=False):
        """ prototype for build-method. Must be implemented in the deriving class """
        raise NotImplementedError("build-method not implemented")

    def __repr__(self):
        """ prototype for _-str__-method. Must be implemented in the deriving class """
        raise NotImplementedError("__rerpr__-method not implemented")


""" ------------------------------------------------------------------------------------
                                        Units
    ------------------------------------------------------------------------------------ """


class Rectangle(ApproximationUnit):
    """ Unit for rectangular-approximation"""

    def __init__(self, x0, width, v):
        """  Initialize the Rectangle-instance

        Parameters
        ----------
        x0 : float
            position of rectangle
        width : float
            width of rectangle
        v : float
            height of the rectangle
        """
        super(Rectangle, self).__init__(x0, width)
        self._v = v

    @property
    def v(self):
        """ float: height of the rectangle """
        return self._v

    def build(self, steps, start=0, zero_edge=False):
        """ Build x- and y-values of the rectangle

        Parameters
        ----------
        steps : int
            x-steps of values
        start : float, optional
            x-offset of values
        zero_edge: bool
            Set edges of rectangle to zero if True

        Returns
        -------
        rectangle_data: tuple of array_like
        """
        x0 = self._x0 + start
        x1 = x0 + self.width
        x = np.linspace(x0, x1, steps)
        y = [self.v] * steps
        if zero_edge:
            y[0] = y[-1] = 0
        return x, y

    def __repr__(self):
        return f"Rectangle(x0={self.x0:.1f}, width={self.width:.1f}, value={self.v:.1f})"
