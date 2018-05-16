from waveprop import Model, Plot, color_gaps
"""
Create model with default leads (v=10, a=1.0, d=0.8).
"""
model = Model()

"""
Creates disordered sample and uses lead-unitcell as prototype.
The disorder stregth is 0.1, so the potential range is v=9-11.
"""
disorder_strength = 0.1  # disorder strength
n_sample_cells = 10  # number of unit-cells in scattering region
model.set_disordered_sample(w=disorder_strength, n=n_sample_cells)

"""
Show the model using included matplotlib-wrapper
"""
plt = model.show()  # plot generated model
plt.show(tight_layout=True)

"""
Plot the transmission curve of the disordered scattering region
"""
energy_range = 0, 30  # energy range to plot the transmission
data = model.transmission_curve(energy_range)  # get transmission data

plt = Plot()

color_gaps(plt, model)  # color band-gaps
plt.plot_transmission(data)  # plot transmission data
plt.show(tight_layout=True)
