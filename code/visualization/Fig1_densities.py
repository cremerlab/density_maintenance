# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import scipy.stats
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv('../../data/literature/collated_total_protein.csv')
mass_spec = pd.read_csv('../../data/literature/collated_mass_fractions_empirics.csv')
delta = 0.0249


# %%
# compute the densities of the mass spec data
membrane = mass_spec[mass_spec['localization'] == 'membrane'].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()
periplasm = mass_spec[mass_spec['localization'] == 'periplasm'].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()

fig, ax = plt.subplots(3, 1, figsize=(2, 2), sharex=True)

# Plot the data
for g, d in size_data.groupby(['source']):
    ax[0].plot(d['growth_rate_hr'], d['volume_um3'], mapper[g]['m'],
               markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=mapper[g]['c'],
               ms=3, alpha=0.75, label=g)
    ax[1].plot(d['growth_rate_hr'], d['surface_to_volume'], mapper[g]['m'],
               markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=mapper[g]['c'],
               ms=3, alpha=0.75)
for g, d in prot_data.groupby(['source']):
    ax[2].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], mapper[g]['m'],
               markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=mapper[g]['c'],
               ms=3, alpha=0.75, label=g)

# # plot the fits
# ax[0].plot(lam_range, vol_fit, lw=1, color=cor['primary_blue'])
# ax[1].plot(lam_range, sav_fit, lw=1, color=cor['primary_blue'])
# ax[2].plot(lam_range, prot_fit, lw=1, color=cor['primary_blue'])

ax[0].set_xlim([0, 2.5])
ax[0].set_ylim([-0.2, 5])
ax[1].set_ylim([3, 10])
ax[2].set_ylim([0, 550])
ax[0].set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
ax[0].set_yticks([0, 1, 3, 5])
ax[2].set_yticks([0, 250, 500])
ax[1].set_yticks([3, 6, 9])
ax[0].set_ylabel('volume\n[µm$^{-3}$]', fontsize=6)
ax[1].set_ylabel('S/V\n[µm$^{-1}$]', fontsize=6)
ax[2].set_ylabel('protein\nper cell [fg]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
# ax[0].legend()
# ax[2].legend()
plt.savefig('../../figures/Fig1_empirical_trends.pdf')

# %%
# Plot the densities
fig, ax = plt.subplots(1, 2, figsize=(3.75, 2), sharex=True)
ax[0].set_ylim([0, 10])
ax[1].set_ylim([0, 175])
for i, _d in enumerate(dfs):
    for g, d in _d.groupby(['dataset_name']):
        ax[i].plot(d['growth_rate_hr'], d['density'], mapper[g]['m'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=mapper[g]['c'],
                   ms=4, alpha=0.5, label=g)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('membrane protein density\n[fg / µm$^2$]', fontsize=6)
ax[1].set_ylabel('periplasm protein density\n[fg / µm$^3$]', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Fig1_empirical_densities.pdf')

# %%
prot_data['volume'] = np.exp(
    vol_popt[1] + vol_popt[0] * prot_data['growth_rate_hr'])
prot_data['density'] = prot_data['fg_protein_per_cell'].values / \
    prot_data['volume']

fig, ax = plt.subplots(1, 1, figsize=(2, 1.5))
for g, d in prot_data.groupby(['source']):
    ax.plot(d['growth_rate_hr'], d['density'], mapper[g]['m'],
            color=mapper[g]['c'], markeredgecolor=cor['primary_black'],
            alpha=0.5, ms=4)

ax.set_ylim([0, 400])
