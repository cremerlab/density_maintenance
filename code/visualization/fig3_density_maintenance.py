#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load our mass spec data
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']

# Load the results of the full inference
samples = pd.read_csv('../analysis/output/phi_rib_scaling_fits_summary.csv')
samples = samples[samples['interval'].isin(['95%', '68%', 'median'])]

# Load the literature mass spec data
lit_data = pd.read_csv('../../data/collated/literature_mass_spec_aggregated.csv')

# Set a mapper for colors
mapper = {'95%':'pale_', '68%':'primary_', 'median':''}
#%%
# Fit linear relations to the data
phi_mem_popt = scipy.stats.linregress(data['phi_rib'], data['phi_mem'])
phi_peri_popt = scipy.stats.linregress(data['phi_rib'], np.log(data['phi_peri']))

phi_rib_range = np.linspace(0, 0.45, 100)
phi_mem_fit = phi_mem_popt.slope * phi_rib_range + phi_mem_popt.intercept
phi_peri_fit = np.exp(phi_peri_popt.slope * phi_rib_range + phi_peri_popt.intercept)

fig, ax = plt.subplots(1, 2, figsize=(3.5, 1.25), sharex=True)
fmt = size.viz.style_point('This Study')


# Plot the literature data
for g, d in lit_data.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.25)
    if g == 'This Study':
        zorder = 1000
    else:
        zorder = 10
    ax[0].plot(d['phi_rib'], d['phi_mem'], zorder=zorder, **fmt)
    ax[1].plot(d['phi_rib'], d['phi_peri'], zorder=zorder, **fmt)

# Plot fits with credible regions
axes = {'phi_mem_ppc': [ax[0], 'blue'], 'phi_peri_ppc': [ax[1], 'purple']}
for g, d in samples[samples['quantity'].isin(['phi_mem_ppc', 'phi_peri_ppc'])].groupby('quantity'):
    _ax, _color = axes[g]
    for _g, _d in d[d['interval'].isin(['95%', '68%', 'median'])].groupby('interval', sort=False):
        if _g == 'median':
            _ax.plot(_d['phi_rib'], _d['lower'], color=cor[f'{mapper[_g]}{_color}'], lw=1)
        else:
            _ax.fill_between(_d['phi_rib'], _d['lower'], _d['upper'], color=cor[f'{mapper[_g]}{_color}'], alpha=0.4)

# Set scaling
ax[0].set_ylim([0, 0.2])
ax[1].set_ylim([0, 0.15])
ax[0].set_xlim([0.05, 0.35])
ax[0].set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax[1].set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax[0].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/fig3_allocation_trends.pdf', bbox_inches='tight')

#%%
# Set kappa values
kappa = [80, 100, 175, 200] 
BETA = 1/0.4558
PHYSICAL_LIMT = 8

fig, ax = plt.subplots(1, 1, figsize=(2, 2))

for k in kappa:
    # Compute theory predictions
    _fit = k * phi_mem_fit / (2 * (1 + BETA * phi_rib_range - phi_mem_fit - phi_peri_fit))

    # Trim to the physical domain
    idx = np.where(_fit <= PHYSICAL_LIMT)[0]
    if len(idx) == 0:
        idx = 0 
    else:
        idx = idx[0]
    ax.plot(phi_rib_range[:idx], _fit[:idx], ':', color=cor['primary_black'], lw=0.5)
    ax.plot(phi_rib_range[idx:], _fit[idx:], '-', color=cor['primary_black'], lw=0.5)

# Plot the predictions
pred = samples[samples['quantity']=='sav_ppc']
for g, d in pred.groupby('interval', sort=False):
    # Do the fill between
    if g == 'median':
        ax.plot(d['phi_rib'], d['lower'], color=cor[f'{mapper[g]}green'], lw=1)
    else:
        ax.fill_between(d['phi_rib'], d['lower'], d['upper'], color=cor[f'{mapper[g]}green'], alpha=0.4)

# Plot the data
fmt = size.viz.style_point('This Study')
ax.plot(data['phi_rib'], data['surface_to_volume'], **fmt)
ax.set_ylim([3.5, 8])
ax.set_xlim([0.05, 0.35])
ax.set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax.set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
ax.set_yticks([4, 5, 6, 7, 8])
ax.set_xticks([0.1, 0.2, 0.3,])
plt.savefig('./plots/fig3_density_maintenance_theory.pdf', bbox_inches='tight')
#%%
# Plot the inferred Kappa
fig, ax = plt.subplots(1, 1, figsize=(2, 1))
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('$\kappa$\ninferred\ndensity ratio [µm$^{-1}$]', fontsize=6)
