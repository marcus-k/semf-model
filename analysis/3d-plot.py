#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# 3D plot of the binding energy per nucleon as a function of N-Z and A.
# 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("../data/binding_energy_per_A.csv")

# Convert from keV to MeV
df["E"] /= 1000        
df["U(E)"] /= 1000

new_df = df.pivot(index="N-Z", columns="A", values="E")

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
xpos, ypos = np.meshgrid(np.arange(new_df.shape[1]), np.arange(new_df.shape[0]))
ax.plot_surface(xpos, ypos, np.exp(new_df.values), cmap='viridis', linewidth=0)

ax.set_xlabel("A")
ax.set_ylabel("N-Z")
ax.set_zlabel("e^E (MeV)")

plt.show()