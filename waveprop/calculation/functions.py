# -*- coding: utf-8 -*-
"""
Created on 13 May 2018
@author: Dylan Jones

general helper methods

"""
import numpy as np
from math import sqrt


def gauss_function(x, a, x0, sigma):
    return a * np.exp(-1/2 * ((x-x0)/sigma)**2)


def mean(array):
    return sum(array)/len(array)


def standard_deviation(array):
    mu = mean(array)
    return sqrt(sum([(x - mu) ** 2 for x in array]) / len(array))


def get_gaps(bands):
    b = bands[0]
    gaps = [[0, b[0]]]
    for i in range(0, len(bands) - 1):
        gaps.append([bands[i][1], bands[i + 1][0]])
    return gaps


def fit_results_mean(results):
    vals, errs = [], []
    for res, err in results:
        vals.append(res)
        errs.append(err)
    return mean(vals), standard_deviation(vals), mean(errs)


def center(val_range):
    return val_range[0] + 0.5*(val_range[1]-val_range[0])


def update_mean(mean_array, array, n_avrg):
    for i in range(len(mean_array)):
        mean_array[i] = mean_array[i] + (array[i] - mean_array[i]) / n_avrg
    return mean_array, n_avrg + 1


def in_region(x, region):
    in_range = (region[0] <= x) and (x <= region[1])
    if in_range:
        return True
