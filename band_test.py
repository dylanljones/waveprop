import numpy as np
import matplotlib.pyplot as plt
from waveprop import build_model_from_data, Cell, Model

lead_cell = Cell.rectangle_barrier(v=8, a=1.0, d=0.8)
sample_cell = Cell.rectangle_barrier(v=10, a=1.2, d=1.0)

model = build_model_from_data()

e_vals = np.linspace(-2, 2, 1000)

for e in e_vals:
  model.set_energy(e)
  lead = model.lead
  print(lead.bloch_vector)