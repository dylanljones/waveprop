# -*- coding: utf-8 -*-
"""
Created on 18 Jul 2018
@author: Dylan Jones

project: wave_propagation
version: 1.0
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from waveprop import build_model_from_data
from waveprop.basic_plots import potential_plot, transmission_plot, transcendental_plot, band_plot

spins = ["up", "down"]
fig_path = r"D:\Dropbox\Work\figures"

save = False
energy_lim = [-30, 30]
energy_steps = 1000


def main():

    for idx in range(3):
        print("----------------------")
        print(f"Loading data-set: {idx:}")
        print("----------------------")
        print("Up: ", end="")
        model_up = build_model_from_data(idx, "up")
        print("Down: ", end="")
        model_down = build_model_from_data(idx, "down")
        models = [model_up, model_down]
        #  e_min = min([min(model.lead.cell.data) for model in models])
        print()
        plot_all(idx, models, energy_lim, energy_steps, save)
        if not save:
            plt.show()
        print()


def plot_all(idx, models, elim, steps=1000, save=False):
    print("  -> Building potential")
    fig1, _ = plot_model(idx, models, 5)
    print("  -> Calculating transcendent equation")
    fig2, _ = plot_transcendent_equation(idx, models, elim, steps)
    print("  -> Calculating band structure")
    fig3, _ = plot_band_structure(idx, models, elim, steps)
    print("  -> Calculating transmission values")
    fig4, _ = plot_transmission(idx, models, elim, steps)
    if save:
        print("Saving Figures")
        dpi = 600
        fig1.savefig(get_fig_name(idx, "a", "model", fig_path), dpi=dpi)
        fig2.savefig(get_fig_name(idx, "b", "transcendental_eq", fig_path), dpi=dpi)
        fig3.savefig(get_fig_name(idx, "c", "bloch_vector", fig_path), dpi=dpi)
        fig4.savefig(get_fig_name(idx, "d", "transmission", fig_path), dpi=dpi)


def get_fig_name(idx, fig_num, name, path=None, ext="png"):
    name = f"{idx:}{fig_num:}_{name:}.{ext:}"
    if path:
        return os.path.join(fig_path, name)
    else:
        return name


def plot_model(idx, models, n):
    v_up, v_down = [model.potential(n) for model in models]
    fig, ax = potential_plot(v_up, v_down)
    ax.set_title("Model (Potential {})".format(idx))
    fig.tight_layout()
    return fig, ax


def plot_transmission(idx, models, elim, steps):
    t_up, t_down = [model.transmission_curve(elim, steps) for model in models]
    fig, ax = transmission_plot(t_up, t_down)
    ax.set_title("Transmission (Potential {})".format(idx))
    ax.set_xlim(*elim)
    fig.tight_layout()
    return fig, ax


def plot_transcendent_equation(idx, models, elim, steps=1000):
    eq_up, eq_down = [model.lead.get_transcendent_equation(elim, steps) for model in models]
    fig, ax = transcendental_plot(eq_up, eq_down)
    ax.set_title("Transcendental Equation (Potential {})".format(idx))
    ax.set_xlim(*elim)
    fig.tight_layout()
    return fig, ax


def plot_band_structure(idx, models, elim, steps=10000):
    band_up, band_down = [model.lead.get_band_structure(elim, steps) for model in models]
    a_up, a_down = [model.lead.cell.a for model in models]
    fig, ax = band_plot(band_up, band_down, a_up, a_down)
    ax.set_title("Bloch Vector (Potential {})".format(idx))
    ax.set_ylim(*elim)
    fig.tight_layout()
    return fig, ax


if __name__ == "__main__":
    main()
