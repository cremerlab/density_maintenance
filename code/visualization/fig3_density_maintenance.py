#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the literature mass spec data
ms_data = pd.read_csv('../../data/collated/literature_mass_spec_aggregated.csv')

# Load our mass spec data
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']

# Load summaries from MCMC
kappa_samples = pd.read_csv('../analysis/output/kappa_density_samples.csv')
pred_sav_summary = pd.read_csv('../analysis/output/predicted_sav_summary.csv')
phi_rib_preds = pd.read_csv('../analysis/output/phi_rib_scaling_fits_summary.csv')
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
physical_maximum = 8
fig, ax = plt.subplots(2, 1, figsize=(2, 4))
fmt = size.viz.style_point('This Study')

for g, d in preds[preds['quantity']=='sav_ppc'].groupby('interval', sort=False):
    idx_lower = np.where(d['lower'] <= 8)[0][0]
    idx_upper = np.where(d['upper'] <= 8)[0][0]
    if g == 'median':
        ax[0].plot(d['phi_rib'].values[:idx_lower], d['lower'][:idx_lower], ':', lw=1, color=cor['red'])
        ax[0].plot(d['phi_rib'].values[idx_lower:], d['lower'][idx_lower:], '-', lw=1, color=cor['green'])
    else:
        ax[0].plot(d['phi_rib'].values[:idx_lower], d['lower'][:idx_lower], ':', lw=1, color=cor[f'{mapper[g]}red'])
        ax[0].plot(d['phi_rib'].values[:idx_upper], d['upper'][:idx_upper], ':', lw=1, color=cor[f'{mapper[g]}red'])

        ax[0].fill_between(d['phi_rib'].values, d['lower'], d['upper'], color=cor[f'{mapper[g]}green'], alpha=0.5)

fmt = size.viz.style_point('This Study')
fmt['alpha'] = 0.95
for g, d in pred_sav_summary[pred_sav_summary['quantity']=='sav_ppc'].groupby(['condition', 'replicate']):
    _d = data[(data['carbon_source']==g[0]) & (data['replicate']==g[1])]    
    for _g, __d in d[d['interval'].isin(['median', '68%', '95%'])].groupby('interval', sort=False):
        if _g == 'median':
            ax[1].plot(_d['surface_to_volume'], __d['lower'],zorder=1000,**fmt)
        else:
            if _g == '68%':
                lw = 1
            else:
                lw = 0.5
            ax[1].vlines(_d['surface_to_volume'], __d['lower'], __d['upper'], color=cor['primary_black'], lw=lw)

ax[1].fill_betweenx([1, 11], 8, 11, color=cor['pale_red'], alpha=0.5)
ax[1].plot([1, 8], [1, 8], '-',color=cor['primary_black'], lw=1)
ax[1].plot([8, 11], [8, 11], ':',color=cor['primary_black'], lw=1)
ax[0].plot(data['phi_rib'], data['surface_to_volume'], **fmt)
ax[0].set_ylim([3, 10])
ax[0].set_xlim([0.03, 0.4])
ax[1].set_xlim([3, 11])
ax[1].set_ylim([3, 11])
ax[0].set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax[0].set_ylabel('$S_A/V$\nsurface-to-volume[µm$^{-1}$]', fontsize=6)
ax[1].set_xlabel('measured surface-to-volume [µm$^{-1}$]', fontsize=6)
ax[1].set_ylabel('predicted surface-to-volume [µm$^{-1}$]', fontsize=6)
plt.savefig('./plots/fig3_theory_test_wildtype.pdf', bbox_inches='tight')

