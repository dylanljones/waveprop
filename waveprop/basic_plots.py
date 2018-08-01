# -*- coding: utf-8 -*-
"""
Created on 1 Aug 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0

This module conatins helper-methods for plotting most relevant curves of model:

    - the potential
    - the transmission
    - the transcendental-equation of the leads
    - the band_structure
"""

import numpy as np
import matplotlib.pyplot as plt


def new_plot(width=16, ratio=5/3):
    """ Creates a new plot with specified size

    Parameters
    ----------
    width: float, optional
        width of the plot in cm
    ratio: float, optional
        ratio of plot

    Returns
    -------
    fig: matplotlib.figure.Figure
    ax:  matplotlib.axes.Axes
    """
    fig, ax = plt.subplots()
    w = width / 2.54
    h = w * ratio
    fig.set_size_inches(h, w)
    return fig, ax


def potential_plot(v_up=None, v_down=None):
    """ Plot potential data

    Parameters
    ----------
    v_up: tuple of array_like, optional
        data of spin-up-potential
    v_down: tuple of array_like, optional
        data of spin-down-potential

    Returns
    -------
    fig: matplotlib.figure.Figure
    ax:  matplotlib.axes.Axes
    """
    fig, ax = new_plot()
    lw = 1.5
    if v_up:
        x, v = v_up
        label = "up"
        ax.plot(x, v, lw=lw, label=label)
    if v_down:
        x, v = v_down
        label = "down"
        ax.plot(x, v, lw=lw, label=label)
    ax.set_xlabel("x [a$_0$]")
    ax.set_ylabel("V [E$_H$]")
    ax.legend()
    fig.tight_layout()
    return fig, ax


def transmission_plot(t_up=None, t_down=None):
    """ Plot the transmission data

    Parameters
    ----------
    t_up: tuple of array_like, optional
        transmission of spin-up-potential
    t_down: tuple of array_like, optional
        transmission of spin-down-potential

    Returns
    -------
    fig: matplotlib.figure.Figure
    ax:  matplotlib.axes.Axes
    """
    fig, ax = new_plot()
    if t_up:
        e, t = t_up
        label = "up"
        ax.plot(e, t, label=label)
    if t_down:
        e, t = t_down
        label = "down"
        ax.plot(e, t, label=label)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("E [E$_H$]")
    ax.set_ylabel("T")
    ax.legend()
    fig.tight_layout()
    return fig, ax


def transcendental_plot(eq_up=None, eq_down=None):
    """ Plot the transcendental equation

    Parameters
    ----------
    eq_up: tuple of array_like, optional
        data of the transcendental equation for the spin-up-potential
    eq_down: tuple of array_like, optional
        data of the transcendental equation for the spin-up-potential

    Returns
    -------
    fig: matplotlib.figure.Figure
    ax:  matplotlib.axes.Axes
    """
    fig, ax = new_plot()
    if eq_up:
        x, y = eq_up
        label = "up"
        ax.plot(x, np.real(y), label=label)
    if eq_down:
        x, y = eq_down
        label = "down"
        ax.plot(x, np.real(y), label=label)
    col = "0.5"
    lw = 0.5
    ls = "--"
    ax.axhline(1, color=col, lw=lw, ls=ls)
    ax.axhline(-1, color=col, lw=lw, ls=ls)

    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("E [E$_H$]")
    ax.set_ylabel("Transcendental Eq.")
    ax.legend()
    fig.tight_layout()
    return fig, ax


def band_plot(band_up=None, band_down=None, a_up=None, a_down=None):
    """ Plot the band-strucutre

    Parameters
    ----------
    band_up: tuple of array_like, optional
        e-, k_re and k_im values of the band_structure of the spin-up-potential
    band_down: tuple of array_like, optional
        e-, k_re and k_im values of the band_structure of the spin-down-potential
    a_up: float, optional
        cell-size of used in the lead of the spin-up-potential
    a_down: float, optional
        cell-size of used in the lead of the spin-up-potential

    Returns
    -------
    fig: matplotlib.figure.Figure
    ax:  matplotlib.axes.Axes
    """
    fig, ax = new_plot()
    i = 0
    if band_up:
        e, k_re, k_im = band_up
        spin = "up"
        col = "C{}".format(i)
        ax.plot(np.real(k_re), e, color=col, label=f"Re({spin:})")
        ax.plot(np.real(k_im), e, color=col, label=f"Im({spin:})", ls="--")
        if a_up:
            ax.axvline(np.pi / a_up, color=col, lw=0.5, ls="--")
            ax.axvline(-np.pi / a_up, color=col, lw=0.5, ls="--")
        i += 1

    if band_down:
        e, k_re, k_im = band_down
        spin = "down"
        col = "C{}".format(i)
        ax.plot(np.real(k_re), e, color=col, label=f"Re({spin:})")
        ax.plot(np.real(k_im), e, color=col, label=f"Im({spin:})", ls="--")
        if a_down:
            ax.axvline(np.pi / a_down, color=col, lw=0.5, ls="--")
            ax.axvline(-np.pi / a_down, color=col, lw=0.5, ls="--")

    ax.set_xlabel("k")
    ax.set_ylabel("E [E$_H$]")
    ax.legend()
    fig.tight_layout()
    return fig, ax
