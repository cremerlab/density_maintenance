# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import size.viz
cor, pal = size.viz.matplotlib_style()
np.random.seed(666)  # Ov reproducibility

# %%
# Generate a heat map of SAV versus width and length
ell_range = np.linspace(0.5, 4, 300)
w_range = np.linspace(0.5, 4, 300)
L, W = np.meshgrid(ell_range, w_range)
SAV = 12 * L / (W * (3 * L - W))
plt.imshow(SAV)
# %%
# Load Literature datasets
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
lit_mass_spec = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')
lit_mass_spec = lit_mass_spec[lit_mass_spec['periplasm'] == True]
lit_mass_spec = lit_mass_spec.groupby(['dataset_name', 'condition',
                                      'growth_rate_hr']
                                      )['mass_frac'].sum().reset_index()

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

prot_data = pd.read_csv(
    '../../data/summaries/summarized_protein_measurements.csv')
prot_data = prot_data[(prot_data['strain'] == 'wildtype') &
                      (prot_data['overexpression'] == 'none') &
                      (prot_data['inducer_conc_ng_mL'] == 0)]
prot_data = prot_data.groupby(['carbon_source'])[
    'mass_frac'].agg(('mean', 'sem')).reset_index()
# Define markers and colors
markers = ['o', 'v', 'X', '<', 's', '>', '^', 'h',
           'p', 'P', '*', 'o', '8', 'd', '>', 'v', '<', '^']
cors = sns.color_palette('Greys_r', n_colors=len(markers)+4).as_hex()[:-4]
np.random.shuffle(cors)

# Get the different data sources and make a mapper
names = list(lit_size_data['source'].unique())
for n in lit_mass_spec['dataset_name'].unique():
    names.append(n)
mapper = {n: {'m': m, 'c': c} for n, m, c in zip(names, markers, cors)}

# %%
# Make the plots of the widths lengths and volumes
fig, ax = plt.subplots(3, 2, figsize=(4, 5))
ax = ax.ravel()
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_ylabel('average length [µm]', fontsize=6)
ax[1].set_ylabel('average width [µm]', fontsize=6)
ax[2].set_ylabel('average volume [µm$^{3}$]', fontsize=6)
ax[3].set_ylabel('surface area to volume [µm$^{-1}$]', fontsize=6)
ax[4].set_ylabel('periplasmic protein\nmass fraction [%]', fontsize=6)
ax[-1].axis('off')
ax[0].set_ylim([1, 6])
ax[1].set_ylim([0.4, 1.3])
ax[2].set_ylim([0, 4.5])
ax[4].set_ylim([0, 15])

for g, d in lit_size_data.groupby(['source']):
    for i, v in enumerate(['length_um', 'width_um', 'volume_um3', 'surface_to_volume']):
        ax[i].plot(d['growth_rate_hr'], d[f'{v}'], linestyle='none', marker=mapper[g]['m'],
                   ms=4, markeredgecolor=cor['primary_black'],
                   markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)
    ax[-1].plot([], [], linestyle='none', marker=mapper[g]['m'],
                ms=4, markeredgecolor=cor['primary_black'],
                markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)


for g, d in lit_mass_spec.groupby(['dataset_name']):
    ax[4].plot(d['growth_rate_hr'], d['mass_frac'] * 100, linestyle='none', marker=mapper[g]['m'],
               ms=4, markeredgecolor=cor['primary_black'],
               markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)
    ax[-1].plot([], [], linestyle='none', marker=mapper[g]['m'],
                ms=4, markeredgecolor=cor['primary_black'],
                markerfacecolor=mapper[g]['c'], alpha=0.5, label=g)


# Plot our data
for g, d in sizes.groupby(['carbon_source']):
    if g == 'ezMOPS':
        continue
    lam = growth_rates[(growth_rates['carbon_source'] == g)]
    for i, v in enumerate(['length', 'width_median', 'volume', 'surface_to_volume']):
        ax[i].errorbar(lam['mean'], d[v]['mean'], xerr=lam['sem'], yerr=d[v]['sem'],
                       marker='o', markeredgecolor=cor['primary_blue'],
                       markerfacecolor='white', markeredgewidth=0.5, ms=4, lw=1,
                       label='__nolegend__', color=cor['primary_blue'])

for g, d in prot_data.groupby(['carbon_source']):
    lam = growth_rates[(growth_rates['carbon_source'] == g)]
    ax[4].errorbar(lam['mean'], d['mean'] * 100, xerr=lam['sem'], yerr=d['sem'], lw=1,
                   marker='o', linestyle='none', markeredgewidth=0.5, ms=4,
                   markerfacecolor='white', markeredgecolor=cor['primary_blue'],
                   label='__nolegend__', color=cor['primary_blue'], capsize=0)

ax[-1].plot([], [], 'o', markeredgecolor=cor['primary_blue'],
            markerfacecolor='white', markeredgewidth=0.5, ms=4, label='This study')
ax[-1].legend()
plt.tight_layout()
plt.savefig('../../figures/Fig1_literature_comparison.pdf',
            bbox_inches='tight')
