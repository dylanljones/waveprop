# -*- coding: utf-8 -*-
"""
Created on 9 Jul 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from ..test_potentials import poeschl_teller_potential
from waveprop.model import Model


def run_poeschl_teller_test():
    elim = 0, 200

    kappas = 1, 1.5, 2, 2.5, 3
    #test_poeschl_teller_lambda(kappas, elim)

    n_rects = 5, 10, 20, 50, 100
    #test_poeschl_teller_rects(n_rects, elim)
    test_poeschl_teller_error(n_rects, elim)


def test_poeschl_teller_lambda(kappas, elim, e_steps=1000):
    fig, ax = plt.subplots()
    for k in kappas:
        e, t = get_transmission(k, elim, n_rects=100, e_steps=e_steps)
        ax.plot(e, t, label="$\lambda$ = {}".format(k))

    ax.legend()
    ax.set_xlim(-0.01*elim[1], elim[1])
    ax.set_ylim(0.01, 1.01)
    ax.set_xlabel("Energy $[E_H]$")
    ax.set_ylabel("Transmission")
    plt.show()


def test_poeschl_teller_rects(n_rects, elim, kappa=1, e_steps=1000):
    fig, ax = plt.subplots()
    for n in n_rects:
        e, t = get_transmission(kappa, elim, n_rects=n, e_steps=e_steps)
        ax.plot(e, t, label="$n_{rect}$ = " + "{}".format(n))
    ax.legend()
    ax.set_xlim(-0.01*elim[1], elim[1])
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlabel("Energy $[E_H]$")
    ax.set_ylabel("Transmission")
    plt.show()


def test_poeschl_teller_error(n_rects, elim, kappa=1):
    fig, ax = plt.subplots()
    for n in n_rects:
        e, t = get_transmission(kappa, elim, n)
        r = 1-np.asarray(t)
        i_mins = argrelextrema(r, np.less)[0]
        e_mins = [e[i] for i in i_mins]
        print(e_mins)

        ax.semilogy(e, r, label="$n_{rect}$ = " + "{}".format(n))

    ax.legend()
    ax.set_xlim(*elim)
    ax.set_xlabel("Energy $[E_H]$")
    ax.set_ylabel("$\ln(1-T)$")
    plt.show()


def build_model(sample_data, lead_data=None, n_rects=20):
    sample_width = abs(sample_data[0][-1] - sample_data[0][0])
    rect_width = sample_width/n_rects
    model = Model()
    model.load_data(sample_data, lead_data)
    model.calculate_rectangles(rect_width)
    return model


def get_transmission(kappa, elim, n_rects, e_steps=1000):
    sample = poeschl_teller_potential(kappa, steps=5000)
    model = build_model(sample, n_rects=n_rects)
    return model.transmission_curve(elim, e_steps)

