# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
# %%
fig, ax = plt.subplots(1, 3, figsize=(6, 3), sharey=True)
y = {'ezMOPS': 6,
     'LB': 5,
     'glucoseCAA': 4,
     'glucose': 3,
     'glycerol': 2,
     'sorbitol': 1,
     'acetate': 0
     }
ax[0].set_xlabel('cell width [µm]')
ax[1].set_xlabel('cell length [µm]')
ax[2].set_xlabel('cell volume [fL]')

ylabs = ['acetate', 'sorbitol', 'glycerol',
         'glucose', 'glucoseCAA', 'LB', 'ezMOPS']
for a in ax.ravel():
    a.set_yticks([0, 1, 2, 3, 4, 5, 6])
    a.set_yticklabels(ylabs)
for g, d in data.groupby(['carbon_source']):
    width = d[d['parameter'] == 'width_um']
    length = d[d['parameter'] == 'length_um']
    volume = d[d['parameter'] == 'volume_fL']
    for i, _d in enumerate([width, length, volume]):
        ax[i].plot([_d['2.5%'], _d['97.5%']], [
            y[g], y[g]], '-', lw=0.5, color='k')
        ax[i].plot([_d['12.5%'], _d['87.5%']], [
            y[g], y[g]], '-', lw=1.5, color='k')
        ax[i].plot([_d['25%'], _d['75%']], [
            y[g], y[g]], '-', lw=3, color='k')
        ax[i].plot(_d['median'], y[g], 'o', ms=2, color='k')
ax[0].set_title('inferred size parameter distributions')
plt.savefig('./wildtype_size_distributions.pdf')
# %%
