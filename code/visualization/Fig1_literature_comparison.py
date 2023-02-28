# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import size.viz
cor, pal = size.viz.matplotlib_style()
np.random.seed(666)  # Ov reproducibility

# Load Literature datasets
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
lit_mass_spec = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')

# Load our datasets
growth_rates = pd.read_csv(
    '../../data/summaries/summarized_growth_measurements.csv')
growth_rates = growth_rates.groupby(['strain', 'carbon_source',
                                     'overexpression', 'inducer_conc']
                                    )['growth_rate_hr'].agg(('mean', 'sem')).reset_index()
growth_rates = growth_rates[(growth_rates['strain'] == 'wildtype') &
                            (growth_rates['overexpression'] == 'none') &
                            (growth_rates['inducer_conc'] == 0)]
sizes = pd.read_csv('../../data/summaries/summarized_size_measurements.csv')
sizes = sizes[(sizes['strain'] == 'wildtype') &
              (sizes['overexpression'] == 'none') &
              (sizes['inducer_conc'] == 0)]
sizes = sizes.groupby(['carbon_source'])[['width_median', 'length',
                                          'volume', 'surface_to_volume']].agg(('mean', 'sem'))

# Define markers and colors
markers = ['o', 'v', 'X', '<', 's', '>', '^', 'h', 'p', 'P', '*', 'o']
cors = sns.color_palette('Greys_r', n_colors=len(markers)+4).as_hex()[:-4]
np.random.shuffle(cors)

# Get the different data sources and make a mapper
names = list(lit_size_data['source'].unique())
for n in lit_mass_spec['dataset_name'].unique():
    names.append(n)
mapper = {n: {'m': m, 'c': c} for n, m, c in zip(names, markers, cors)}

# %%
# Make the plots of the widths lengths and volumes
fig, ax = plt.subplots(1, 3, figsize=(5, 1.5))
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_ylabel('average length [µm]', fontsize=6)
ax[1].set_ylabel('average width [µm]', fontsize=6)
ax[2].set_ylabel('average volume [µm$^{3}$]', fontsize=6)
ax[0].set_ylim([0, 6])
ax[1].set_ylim([0, 1.5])
ax[2].set_ylim([0, 4.5])

for g, d in lit_size_data.groupby(['source']):
    for i, v in enumerate(['length_um', 'width_um', 'volume_um3']):
        ax[i].plot(d['growth_rate_hr'], d[f'{v}'], linestyle='none', marker=mapper[g]['m'],
                   ms=4, markeredgecolor=cor['primary_black'],
                   markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)

# Plot our data
for g, d in sizes.groupby(['carbon_source']):
    if g == 'ezMOPS':
        continue
    lam = growth_rates[(growth_rates['carbon_source'] == g)]
    for i, v in enumerate(['length', 'width_median', 'volume']):
        ax[i].errorbar(lam['mean'], d[v]['mean'], xerr=lam['sem'], yerr=d[v]['sem'], marker='D', markeredgecolor=cor['primary_green'],
                       markerfacecolor='white', markeredgewidth=1, ms=4, lw=1, label='__nolegend__')

ax[0].plot([], [], 'D', markeredgecolor=cor['primary_green'],
           markerfacecolor='white', markeredgewidth=1, ms=4, label='This study')
ax[0].legend()
plt.tight_layout()
plt.savefig('../../figures/Fig1_literature_comparison.pdf',
            bbox_inches='tight')

# %%
fig, ax = plt.subplots(2, 2,)
ax = ax.ravel()
ax[1].set_ylim([1, 6])

for g, d in lit_size_data.groupby(['source']):
    # for i, v in enumerate(['length_um', 'width_um', 'volume_um3']):
    ax[0].plot(d['growth_rate_hr'], d['surface_to_volume'], linestyle='none', marker=mapper[g]['m'],
               ms=4, markeredgecolor=cor['primary_black'],
               markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)
    ax[1].plot(d['growth_rate_hr'], d['length_um'] / d['width_um'], linestyle='none', marker=mapper[g]['m'],
               ms=4, markeredgecolor=cor['primary_black'],
               markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)
