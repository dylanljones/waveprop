# -*- coding: utf-8 -*-
"""
Created on 9 May 2018
@author: Dylan Jones

This module contains the Constants object, wich gets initialized for later use

"""
from scipy import constants as c

_m = 1  # mass of the particle: Set to one for atomic units
_hbar = 1  # hbar: Set to one for atomic units


class Constants:
    """ Constants object containing all used constants"""

    def __init__(self):
        self.m = _m
        self.hbar = _hbar
        self.e_h = c.hbar * c.c * c.alpha * c.physical_constants["Bohr radius"][0]


constants = Constants()  # create the constants object for later use
