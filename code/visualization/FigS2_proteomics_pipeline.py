# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

# Load the data sets
ms_raw = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ms_emp = pd.read_csv('../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
lam_fits = pd.read_csv(
    '../../data/mcmc/theory_growth_rate_prediction_summaries.csv')


# Data Panel 1: Raw MS fractions
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.set_xlim([0, 2.5])
ax.set_ylim([0, 0.25])
ax.set_ylabel('$\phi_{mem}$\nmembrane proteome fraction', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

for g, d in ms_raw[ms_raw['localization'] == 'membrane'].groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['light_blue']
    ax.plot(d['growth_rate_hr'], d['mass_frac'], **fmt, ms=4)
ax.legend(fontsize=6)
plt.savefig('../../figures/FigS2_membrane_fraction.pdf')


# %%

# Data Panel 2: Surface area and protein per cell
fig, ax = plt.subplots(2, 1, figsize=(1.5, 1.5), sharex=True)
ax[1].set_xticks([0, 0.5, 1, 1.5, 2, 2.5])
ax[1].set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[0].set_ylim([0, 800])
ax[0].set_yticks([0, 200, 400, 600, 800])
ax[0].set_ylabel('$M_{prot}^{(tot)}$ [fg / cell]\ntotal protein', fontsize=6)
for g, d in prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt, ms=4)

ax[1].set_ylim([1, 10])
ax[1].set_yticks([2, 4, 6, 8, 10])
ax[1].set_ylabel('$S_A$ [µm$^2$]\nsurface area', fontsize=6)
for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt, ms=4)

axes = {'pred_lam_prot': ax[0], 'SA_fit': ax[1]}
inter_colors = {'95%': 'pale_', '75%': 'light_',
                '25%': 'primary_', 'median': ''}
for g, d in lam_fits[lam_fits['quantity'].isin(list(axes.keys()))].groupby(['quantity', 'interval'], sort=False):
    axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], alpha=0.5,
                            color=cor[f'{inter_colors[g[1]]}blue'])

ax[0].legend(fontsize=6)
ax[1].legend(fontsize=6)
plt.savefig('../../figures/FigS2_prot_sa_fits.pdf')

# %%

# Data Panel 1: Raw MS fractions
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.set_xlim([0, 2.5])
ax.set_ylim([0, 6])
ax.set_ylabel(r'$\rho_{mem}$ [fg / µm$^{2}$]' +
              '\nmembrane protein density', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

for g, d in ms_emp[ms_emp['quantity'] == 'ms_rho_mem'].groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['light_blue']
    ax.vlines(d['growth_rate_hr'], d['2.5%'],
              d['97.5%'], color=fmt['color'], lw=0.5)
    ax.plot(d['growth_rate_hr'], d['median_value'], **fmt, ms=4, zorder=1000)

plt.savefig('../../figures/FigS2_membrane_density.pdf')
