# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()
ppc_cmap = {'95%': cor['pale_black'], '75%': cor['light_black'],
            '50%': cor['primary_black'], '25%': cor['black'], '10%': cor['dark_black']}

# Load the various datasets
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data['aspect_ratio'] = size_data['length_um'].values / \
    size_data['width_um'].values
ppcs = pd.read_csv('../../data/mcmc/literature_model_params.csv')
ppcs = ppcs[ppcs['model'] == 'const_phi_mem']


# %%

fig, ax = plt.subplots(3, 1, figsize=(2, 2.25), sharex=True)
ax[0].set_ylim([0, 20])
ax[1].set_ylim([0, 0.25])
ax[2].set_ylim([1, 8])

# Isolate the median values for the constants
const_medians = ppcs[(ppcs['interval'] == 'median') &
                     (ppcs['quantity'].isin(['alpha', 'm_peri_mu', 'phi_mem_mu']))]

# Plot the literature data
for g, d in ms_data[ms_data['localization'] == 'periplasm'].groupby(['dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['mass_fg'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
               markeredgewidth=0.5, alpha=0.5)

for g, d in ms_data[ms_data['localization'] == 'membrane'].groupby(['dataset_name']):
    ax[1].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
               markeredgewidth=0.5, alpha=0.5)

for g, d in size_data.groupby(['source']):
    ax[2].plot(d['growth_rate_hr'], d['aspect_ratio'], mapper[g]['m'],
               ms=3, markeredgecolor=cor['primary_black'], color=mapper[g]['c'],
               markeredgewidth=0.5, alpha=0.5)

axes = {ax[0]: 'm_peri_mu', ax[1]: 'phi_mem_mu', ax[2]: 'alpha_mu'}
for g, d in const_medians.groupby(['quantity'], sort=False):
    d.sort_values(by='growth_rate_hr')
    axes[g].plot(d['growth_rate_hr'], d['lower'], '-',
                 lw=1.5, color=cor['primary_blue'])
