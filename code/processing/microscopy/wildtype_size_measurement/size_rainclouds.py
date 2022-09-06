# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

# Load all of the cell sizes
cell_sizes = pd.read_csv('./output/wildtype_size_measurements.csv')
cell_sizes
# %%
cell_sizes['carbon_source'].unique()
# %%
fig, ax = plt.subplots(
    len(cell_sizes['carbon_source'].unique()), 1, figsize=(6, 6), sharex=True)
i = 0
for g, d in cell_sizes.groupby(['carbon_source']):
    for _g, _d in d.groupby(['date']):
        ax[i].hist(_d['width_median'], bins=20,
                   label=g, alpha=0.75, density=True)
    ax[i].set_title(g)
    i += 1

# %%
