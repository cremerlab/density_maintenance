# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import size.viz
cor, pal = size.viz.matplotlib_style()
np.random.seed(666)  # Ov reproducibility

# Load the datasets
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
mass_spec = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
params = params[(params['strain'] == 'wildtype') & (params['overexpression'] == 'none') &
                (params['inducer_conc'] == 0)]
param_medians = params[params['interval'] == 'median']
param_errs = params[params['interval'] != 'median']

# Define markers and colors
markers = ['o', 'v', 'X', '<', 's', '>', '^', 'h', 'p', 'P', '*', 'o']
cors = sns.color_palette('Greys_r', n_colors=len(markers)+4).as_hex()[:-4]
np.random.shuffle(cors)

# Get the different data sources and make a mapper
names = list(size_data['source'].unique())
for n in mass_spec['dataset_name'].unique():
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
plt.tight_layout()

for g, d in size_data.groupby(['source']):
    for i, v in enumerate(['length_um', 'width_um', 'volume_um3']):
        ax[i].plot(d['growth_rate_hr'], d[f'{v}'], linestyle='none', marker=mapper[g]['m'],
                   ms=4, markeredgecolor=cor['primary_black'],
                   markerfacecolor=mapper[g]['c'], alpha=0.5)

# Plot our data
err_lw = {'95%': 0.5, '75%': 1,  '25%': 1.5}
for g, d in param_medians.groupby(['carbon_source']):
    lam = d[d['quantity'] == 'growth_mu']
    for i, v in enumerate(['length_mu', 'width_mu', 'volume_mu']):
        ax[i].plot(lam['lower'], d[d['quantity'] == v]['lower'], 'D', markeredgecolor=cor['green'],
                   markerfacecolor=cor['pale_green'], markeredgewidth=1, ms=4)
