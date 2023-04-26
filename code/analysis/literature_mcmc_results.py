# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()
ppc_cmap = {'95%': cor['pale_blue'], '75%': cor['light_blue'],
            '50%': cor['primary_blue'], '25%': cor['blue'],
            '10%': cor['dark_blue']}

# Load the various datasets
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv(
    '../../data/literature/collated_total_protein_density.csv')

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ms_data['surface_area'] = ms_data['surface_to_volume'] * ms_data['volume']

# Load the various ppcs
ppcs = pd.read_csv('../../data/mcmc/literature_model_params.csv')
ppcs
mode = 'rep'
model = 'const_phi_mem'
ppcs = ppcs[ppcs['model'] == model]

fig, ax = plt.subplots(1, 4, figsize=(6, 1.5))

# Plot the posterior predictive checks
ax = ax.ravel()
axes = {'w': ax[0], 'ell': ax[1], 'vol': ax[2]}
for g, d in ppcs[(ppcs['interval'] != 'median') & (ppcs['quantity'].isin([f'{p}_{mode}' for p in ['w', 'ell', 'vol']]))].groupby(['quantity', 'interval'], sort=False):
    d.sort_values(by='growth_rate_hr', inplace=True)
    axes[g[0].split('_')[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=ppc_cmap[g[1]],
                                          alpha=0.3)

for g, d in ppcs[(ppcs['interval'] != 'median') & (ppcs['quantity'] == f'prot_per_cell_{mode}')].groupby(['quantity', 'interval'], sort=False):
    d.sort_values(by='growth_rate_hr', inplace=True)
    ax[3].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=ppc_cmap[g[1]],
                       alpha=0.3)

# Plot the data
for g, d in size_data.groupby(['source']):
    for i, v in enumerate(['width_um', 'length_um', 'volume_um3']):
        ax[i].plot(d['growth_rate_hr'], d[f'{v}'], mapper[g]['m'], color=mapper[g]['c'],
                   alpha=0.5, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   ms=4)

# Plot the total protein data
for g, d in prot_data.groupby(['source']):
    ax[3].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], mapper[g]['m'], color=mapper[g]['c'],
               alpha=0.5, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
               ms=4)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_xlim([0, 2.5])
ax[2].set_ylim([0, 5])
ax[3].set_ylim([0, 500])
ax[3].set_xlim([0, 2.5])

ax[0].set_ylabel('average width [µm]', fontsize=6)
ax[1].set_ylabel('average length [µm]', fontsize=6)
ax[2].set_ylabel('average volume [µm$^3$]', fontsize=6)
ax[3].set_ylabel('protein per cell [fg]', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Fig2_size_prot_ppcs.pdf', bbox_inches='tight')

# %%
total_prot = ms_data[ms_data['localization'].isin(['cytoplasm',
                                                   'envelope'])
                     ].groupby(['dataset_name', 'condition', 'growth_rate_hr',
                                'volume'])['mass_fg'].sum().reset_index()

# Plot the allocation ppcs
fig, ax = plt.subplots(3, 2, figsize=(6, 4.5))
ax = ax.ravel()
axes = {'membrane': [ax[0], ax[2]], 'periplasm': [ax[1], ax[3]]}
for g, d in ms_data[ms_data['localization'].isin(['membrane', 'periplasm', 'cytoplasm'])].groupby(['dataset_name', 'localization']):
    if g[1] == 'membrane':
        d['density'] = d['mass_fg'] / (2 * d['surface_area'])
    elif g[1] == 'periplasm':
        d['density'] = d['mass_fg'] / (d['surface_area'] * 0.0249)
    else:
        d = total_prot[total_prot['dataset_name'] == g[0]]
        d['density'] = d['mass_fg'] / d['volume']
    if g[1] == 'cytoplasm':
        ax[5].plot(d['growth_rate_hr'], d['density'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
                   alpha=0.5, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   ms=4)
    else:
        axes[g[1]][0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
                           alpha=0.5, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                           ms=4)
        axes[g[1]][1].plot(d['growth_rate_hr'], d['density'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
                           alpha=0.5, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                           ms=4)

        if g[1] == 'periplasm':
            ax[-2].plot(d['growth_rate_hr'], d['mass_fg'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
                        alpha=0.5, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                        ms=4)

axes = {'phi_mem': ax[0], 'phi_peri': ax[1],
        'rho_mem': ax[2], 'rho_peri': ax[3], 'm_peri': ax[4],
        'rho_prot': ax[5]}
for g, d in ppcs[(ppcs['interval'] != 'median') &
                 (ppcs['quantity'].isin([f'{p}_{mode}' for p in axes.keys()]))].groupby(['quantity', 'interval'], sort=False):
    d.sort_values('growth_rate_hr', inplace=True)
    axes[g[0].split(f'_{mode}')[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                                                 color=ppc_cmap[g[1]], alpha=0.3)

ax[0].set_ylim([0, 0.2])
ax[1].set_ylim([0, 0.125])
ax[2].set_ylim([0, 6])
ax[3].set_ylim([0, 175])
ax[4].set_ylim([0, 20])
for a in ax:
    a.set_xlim([0, 2])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('allocation toward\nmembrane protein', fontsize=6)
ax[1].set_ylabel('allocation toward\nperiplasm protein', fontsize=6)
ax[2].set_ylabel('membrane protein density\n[fg / µm$^2$]', fontsize=6)
ax[3].set_ylabel('periplasmic protein density\n[fg / µm$^3$]', fontsize=6)
ax[4].set_ylabel('periplasmic protein mass [fg]', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Fig2_model_trends.pdf', bbox_inches='tight')

# %%
ppcs = pd.read_csv('../../data/mcmc/literature_model_params.csv')
fig, ax = plt.subplots(1, 2, figsize=(4, 2), sharey=True)
ax[0].set_xlim([0.5, 1.1])
for a in ax:
    a.set_ylim([0, 1])
    a.set_xlabel('average cell width [µm]', fontsize=6)
ax[0].set_ylabel(
    'relative allocation of envelope\n$\phi_{peri}/\phi_{mem}$', fontsize=6)
for g, d in ms_data[ms_data['localization'].isin(['membrane', 'periplasm'])].groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    rel_alloc = d[d['localization'] == 'periplasm']['mass_frac'].values[0] / \
        d[d['localization'] == 'membrane']['mass_frac'].values[0]
    for i in range(2):
        ax[i].plot(d['width'].values[0], rel_alloc,
                   mapper[g[0]]['m'], color=mapper[g[0]]['c'],
                   markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   ms=4, alpha=0.5)

const_rho = ppcs[ppcs['model'] == 'const_rho_mem']
const_phi = ppcs[ppcs['model'] == 'const_phi_mem']
for g, d in const_phi[(const_phi['interval'] != 'median') & (const_phi['quantity'].isin([f'rel_phi_{mode}']))].groupby(['interval'], sort=False):
    d.sort_values(by='width', inplace=True)
    ax[0].fill_between(d['width'], d['lower'], d['upper'], color=ppc_cmap[g],
                       alpha=0.3)

for g, d in const_rho[(const_rho['interval'] != 'median') & (const_rho['quantity'].isin([f'rel_phi_{mode}']))].groupby(['interval'], sort=False):
    d.sort_values(by='width', inplace=True)
    ax[1].fill_between(d['width'], d['lower'], d['upper'], color=ppc_cmap[g],
                       alpha=0.3)

ax[1].set_xlim([0.5, 1.25])
