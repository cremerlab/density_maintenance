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
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')

# Load the various ppcs
size_ppcs = pd.read_csv('../../data/mcmc/literature_model_size_ppcs.csv')
prot_ppcs = pd.read_csv('../../data/mcmc/literature_model_protein_ppcs.csv')
ms_ppcs = pd.read_csv('../../data/mcmc/literature_model_ms_ppcs.csv')

# %%
fig, ax = plt.subplots(1, 4, figsize=(6, 1.5))

# Plot the posterior predictive checks
ax = ax.ravel()
axes = {'w_rep': ax[0], 'ell_rep': ax[1], 'vol_rep': ax[2]}
for g, d in size_ppcs[size_ppcs['interval'] != 'median'].groupby(['quantity', 'interval'], sort=False):
    d.sort_values(by='growth_rate_hr', inplace=True)
    axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=ppc_cmap[g[1]],
                            alpha=0.5)

for g, d in prot_ppcs[prot_ppcs['interval'] != 'median'].groupby(['quantity', 'interval'], sort=False):
    d.sort_values(by='growth_rate_hr', inplace=True)
    ax[3].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=ppc_cmap[g[1]],
                       alpha=0.5)

# Plot the data
for g, d in size_data.groupby(['source']):
    for i, v in enumerate(['width_um', 'length_um', 'volume_um3']):
        ax[i].plot(d['growth_rate_hr'], d[f'{v}'], mapper[g]['m'], color=mapper[g]['c'],
                   alpha=0.75, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                   ms=4)

# Plot the total protein data
for g, d in prot_data.groupby(['source']):
    ax[3].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], mapper[g]['m'], color=mapper[g]['c'],
               alpha=0.75, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
               ms=4)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_xlim([0, 2.5])
ax[2].set_ylim([0, 5])
ax[3].set_ylim([0, 550])

ax[0].set_ylabel('average width [µm]', fontsize=6)
ax[1].set_ylabel('average length [µm]', fontsize=6)
ax[2].set_ylabel('average volume [µm$^3$]', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Fig2_size_prot_ppcs.pdf', bbox_inches='tight')

# %%
# Plot the allocation ppcs
fig, ax = plt.subplots(2, 2, figsize=(6, 3))
ax = ax.ravel()
axes = {'membrane': [ax[0], ax[2]], 'periplasm': [ax[1], ax[3]]}
for g, d in ms_data[ms_data['localization'].isin(['membrane', 'periplasm'])].groupby(['dataset_name', 'localization']):
    axes[g[1]][0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g[0]]['m'], color=mapper[g[0]]['c'],
                       alpha=0.75, markeredgecolor=cor['primary_black'], markeredgewidth=0.5,
                       ms=4)

# # Plot the posterior predictive checks
# ax = ax.ravel()
# axes = {'w_rep': ax[0], 'ell_rep': ax[1], 'vol_rep': ax[2]}
# for g, d in size_ppcs[size_ppcs['interval'] != 'median'].groupby(['quantity', 'interval'], sort=False):
#     d.sort_values(by='growth_rate_hr', inplace=True)
#     axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=ppc_cmap[g[1]],
#                             alpha=0.5)