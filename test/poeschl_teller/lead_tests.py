# -*- coding: utf-8 -*-
"""
Created on 15 Jul 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0
"""
import matplotlib.pyplot as plt
from waveprop import Cell, Model, Approximation
from ..test_potentials import poeschl_teller_potential


def run_poeschl_teller_lead():
    sample_width = 4
    kappa = 1
    sample_data = poeschl_teller_potential(kappa+1, sample_width, sample_width)
    sample_approx = Approximation(*sample_data)
    sample_approx.rectangles(0.2)

    cell_data = poeschl_teller_potential(kappa, sample_width, sample_width)
    cell_approx = Approximation(*cell_data)
    cell_approx.rectangles(0.2)
    cell = Cell(cell_approx.units)

    cell = Cell.rectangle_barrier(2, 3, 1)

    model = Model()
    model.load_approximation(sample_approx)
    model.load_lead(cell)

    fig1 = plot_model(model)
    fig2 = plot_transmission(model)

    plt.show()


def plot_model(model):
    fig, ax = new_plot(16)
    ax.set_title("Pöschl-Teller Lead")
    x, v = model.potential(4)
    ax.plot(x, v)
    ax.set_xlabel("x [a$_0$]")
    ax.set_ylabel("V [E$_H$]")
    fig.tight_layout()
    return fig


def plot_transmission(model):
    fig, ax = new_plot(16)
    ax.set_title("Pöschl-Teller Lead")
    e, t = model.transmission_curve([0, 20])
    ax.plot(e, t)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("E [E$_H$]")
    ax.set_ylabel("T")
    fig.tight_layout()
    return fig


def new_plot(width=12, ratio=6/4):
    fig, ax = plt.subplots()
    w = width / 2.54
    h = w * ratio
    fig.set_size_inches(h, w)
    return fig, ax

