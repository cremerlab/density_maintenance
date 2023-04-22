# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv(
    '../../data/literature/collated_total_protein_density.csv')
mass_spec = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
fit = pd.read_csv('../../data/empirical_literature_trends.csv')
delta = 0.0249
# prot_data = prot_data[prot_data['source'] != 'Dennis & Bremer 1974']
# prot_data = prot_data[prot_data['source'] != 'Chohji et al. 1976']

# %%
# compute the densities of the mass spec data
membrane = mass_spec[mass_spec['localization'] == 'membrane'].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'surface_to_volume', 'volume'])['mass_fg'].sum().reset_index()
periplasm = mass_spec[mass_spec['localization'] == 'periplasm'].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'surface_to_volume', 'volume'])['mass_fg'].sum().reset_index()
cyto = mass_spec[mass_spec['localization'] == 'cytoplasm'].groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'surface_to_volume', 'volume'])['mass_fg'].sum().reset_index()

fig, ax = plt.subplots(3, 1, figsize=(2, 2.5), sharex=True)

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
ax[0].plot(fit['growth_rate_hr'], fit['volume'],
           lw=1, color=cor['primary_blue'])
ax[1].plot(fit['growth_rate_hr'], fit['surface_to_volume'],
           lw=1, color=cor['primary_blue'])
ax[2].plot(fit['growth_rate_hr'], fit['fg_protein_per_cell'],
           lw=1, color=cor['primary_blue'])
ax[0].legend()
ax[2].legend()
ax[0].set_xlim([0, 2.5])
ax[0].set_ylim([-0.2, 5])
ax[1].set_ylim([3, 10])
ax[2].set_ylim([0, 800])
ax[0].set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
ax[0].set_yticks([0, 1, 3, 5])
# ax[2].set_yticks([0, 250, 500])
ax[1].set_yticks([3, 6, 9])
ax[0].set_ylabel('volume\n[µm$^{-3}$]', fontsize=6)
ax[1].set_ylabel('S/V\n[µm$^{-1}$]', fontsize=6)
ax[2].set_ylabel('protein\nper cell [fg]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
plt.savefig('../../figures/Fig1_empirical_trends.pdf')
# %%
# Plot the densities
fig, ax = plt.subplots(2, 2, figsize=(3.75, 2.5), sharex=True)
ax = ax.ravel()
ax[0].set_ylim([150, 250])
ax[1].set_ylim([0, 6])
ax[2].set_ylim([0, 170])
ax[3].set_ylim([0.5, 7])
ax[3].set_yticks([1, 3, 5, 7])
density_colors = [cor['primary_black'],
                  cor['light_blue'], cor['primary_purple']]
for i, _d in enumerate([cyto, membrane, periplasm]):
    if i == 1:
        _d['density'] = _d['mass_fg'] / \
            (2 * _d['surface_to_volume'] * _d['volume'])
    elif i == 2:
        _d['density'] = _d['mass_fg'] / \
            (_d['surface_to_volume'] * _d['volume'] * delta)
    else:
        _d['density'] = _d['mass_fg'] / _d['volume']
    _d.sort_values(by='growth_rate_hr', inplace=True)
    _d['fc'] = _d['density'] / _d['density'].min()
    for g, d in _d.groupby(['dataset_name']):

        ax[i].plot(d['growth_rate_hr'], d['density'], mapper[g]['m'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=density_colors[i],
                   ms=4, alpha=0.5, label=g)

        ax[3].plot(d['growth_rate_hr'], d['fc'], mapper[g]['m'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=density_colors[i],
                   ms=4, alpha=0.5, label=g)

for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_ylabel(r'$\rho_{cyt}$'+'\ndensity [fg / µm$^3$]', fontsize=6)
ax[0].set_title('cytoplasmic protein mass', fontsize=6)
ax[1].set_ylabel(r'$\rho_{mem}$'+'\nareal density [fg / µm$^2$]', fontsize=6)
ax[1].set_title('membrane protein mass', fontsize=6)
ax[2].set_ylabel(r'$\rho_{peri}$'+'\ndensity [fg / µm$^3$]', fontsize=6)
ax[2].set_title('periplasmic protein mass', fontsize=6)
ax[3].set_ylabel(
    r'$\rho / \rho_{min}$', fontsize=6)
ax[3].set_title('fold-change in density', fontsize=6)
# plt.tight_layout()
plt.savefig('../../figures/Fig1_empirical_densities.pdf')
# %%

# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_ylim([0, 500])
for g, d in prot_data.groupby(['source']):
    plt.plot(d['growth_rate_hr'], d['density'], mapper[g]['m'],
             markeredgecolor=cor['primary_black'], markeredgewidth=0.5, color=mapper[g]['c'],
             ms=3, alpha=0.75, label=g)

# ax.legend()
