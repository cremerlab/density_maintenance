# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Load the mass spec data
data = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
data = data[data['dataset_name'] != 'Valgepea et al. 2013']

cyto = data[data['envelope'] == False]
cyto = cyto.groupby(['dataset_name', 'condition',
                     'growth_rate_hr', 'ribosomal', 'metabolism'])['mass_frac'].sum().reset_index()
membrane = data[data['membrane'] == True]
membrane = membrane.groupby(['dataset_name', 'condition',
                             'growth_rate_hr'])['mass_frac'].sum().reset_index()
periplasm = data[data['periplasm'] == True]
periplasm = periplasm.groupby(['dataset_name', 'condition',
                               'growth_rate_hr'])['mass_frac'].sum().reset_index()

fig, ax = plt.subplots(1, 3, figsize=(6, 1.5))
for g, d in membrane.groupby(['dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=cor['dark_purple'], alpha=0.5)

for g, d in periplasm.groupby(['dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'], ms=4,
               markeredgecolor='k', markeredgewidth=0.5, color=cor['light_purple'], alpha=0.5)

for g, d in data[data['envelope'] == True].groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    peri = d[d['periplasm'] == True]
    mem = d[d['membrane'] == True]
    ax[2].plot(g[-1], peri['mass_frac'].sum() / mem['mass_frac'].sum(), mapper[g[0]]['m'],
               markerfacecolor=cor['purple'], markeredgecolor=cor['primary_black'],
               markeredgewidth=0.5, ms=4, alpha=0.75)
ax[0].set_ylim([0, 0.25])
ax[1].set_ylim([0, 0.15])
ax[2].set_ylim([0, 1])

# %%

# %%

# %%
