# -*- coding: utf-8 -*-
"""
Created on 29 Jun 2018
@author: Dylan Jones

project: first_try
version: 1.0
"""
import numpy as np


def round_potential(a, v, v0=0, steps=200):
    x = np.linspace(0, a, steps)
    y = 0.5 * v * (np.cos(2 * np.pi / a * (x - a / 2)) + 1) + v0
    return x, y


def poeschl_teller_potential(kappa, x0=0, x_range=4, steps=200):
    x = np.linspace(-x_range, x_range, steps)
    y = -0.5 * (kappa * (kappa + 1)) * (np.cosh(x))**(-2)
    return x + x0, y
