#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

ms_data = pd.read_csv('../../data/collated/literature_mass_spec_aggregated.csv')
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']
preds = pd.read_csv('../analysis/output/phi_rib_scaling_fits_summary.csv')
lit_data = pd.read_csv('../../data/literature/Si2017/si2017_size_phi_rib.csv')
intervals = ['95%', '68%', 'median']
preds = preds[preds['interval'].isin(intervals)]

#%%
fig, ax = plt.subplots(2, 1, figsize=(2, 2), sharex=True)

# Plot the percentiles
quantities = ['phi_mem_ppc', 'phi_peri_ppc']
mapper = {'95%':'pale_', '68%':'primary_', 'median':''}
for i, q in enumerate(quantities):
    for g, d in preds[preds['quantity']==q].groupby('interval', sort=False):
        if i == 0:
            color = cor[f'{mapper[g]}blue']
        else:
            color = cor[f'{mapper[g]}purple']

        if g != 'median':
            ax[i].fill_between(d['phi_rib'], d['lower'], d['upper'], color=color, alpha=0.5)
        else:
            ax[i].plot(d['phi_rib'], d['lower'], lw=1, color=color)

for g, d in ms_data.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.25)
    ax[0].plot(d['phi_rib'], d['phi_mem'], **fmt)
    ax[1].plot(d['phi_rib'], d['phi_peri'], **fmt)

# Set the axis limits and labels
ax[0].set_ylim([0.05, 0.20])
ax[0].set_xlim([0.05, 0.35])
ax[1].set_ylim([0, 0.15])
ax[1].set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax[0].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)

#%%


#%%
fig, ax = plt.subplots(2, 1, figsize=(2, 3))
si_fmt = size.viz.style_point('Si et al. 2017')
fmt = size.viz.style_point('This Study')
for g, d in preds[preds['quantity']=='sav_ppc'].groupby('interval', sort=False):
    if g == 'median':
        ax[0].plot(d['phi_rib'], d['lower'], lw=1, color=cor['green'])
    else:
        ax[0].fill_between(d['phi_rib'], d['lower'], d['upper'], color=cor[f'{mapper[g]}green'], alpha=0.2)
ax[0].plot(lit_data['phi_rib'], lit_data['surface_to_volume'],**si_fmt)
ax[0].plot(data['phi_rib'], data['surface_to_volume'], **fmt)
ax[0].set_ylim([2, 9])

