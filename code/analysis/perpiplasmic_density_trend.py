# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import size.analytical
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/mcmc/wildtype_hyperparameter_size_samples.csv')
voldata = pd.read_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
grdata = pd.read_csv('../../data/mcmc/wildtype_growth_rate_summary.csv')
# %%
_dfs = []
data['periplasmic_volume_fL'] = size.analytical.envelope_volume(
    data['length_um'], data['width_um'], 0.025)

# Define the percentiles to compute for each hyperparameter
percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']

growth_rates = {'acetate': {'growth_rate': 0.35,
                            'protein': 23},
                'sorbitol': {'growth_rate': 0.46,
                             'protein': 18},
                'glycerol': {'growth_rate': 0.68,
                             'protein': 16.8},
                'glucose': {'growth_rate': 0.87,
                            'protein': 10.76},
                'glucoseCAA': {'growth_rate': 1.17,
                               'protein': 10.33}}


def cell_number(vol):
    return 1/(1.3372E-9 * vol - 4.71E-10)


for g, d in data.groupby(['carbon_source']):
    if g in growth_rates.keys():
        # Compute various summary statistics of the hyper parameters
        prot_per_cell = growth_rates[g]['protein'] / \
            cell_number(d['volume_fL']) * 1E9
        density = prot_per_cell / d['periplasmic_volume_fL']
        _percs = np.percentile(density, percs)
        _perc_df = pd.DataFrame([_percs], columns=perc_cols)
        _perc_df['mean'] = density.mean()
        _perc_df['median'] = density.median()
        _perc_df['carbon_source'] = g
        _perc_df['strain'] = 'wildtype'
        _dfs.append(_perc_df)
df = pd.concat(_dfs, sort=False)
# %%

mass_spec = pd.read_csv(
    '../../data/calculated_envelope_protein_mass_and_concentration.csv')
fig, ax = plt.subplots(1, 2, figsize=(6, 4))
# plt.plot(mass_spec['growth_rate_hr'],
#  mass_spec['periplasmic_protein_density_fg_fL'], '.')
# plt.plot(df['growth_rate'], d['mean'], 'ko')
for g, d in df.groupby(['carbon_source']):
    _vol = voldata[(voldata['carbon_source'] == g) &
                   (voldata['parameter'] == 'volume_fL')]
    _gr = grdata[grdata['carbon_source'] == g]
    ax[1].vlines(_gr['median'], d['2.5%'],
                 d['97.5%'], lw=1, color='k')
    ax[1].vlines(_gr['median'], d['12.5%'],
                 d['87.5%'], lw=2, color='k')
    ax[1].hlines(d['median'], _gr['2.5%'],
                 _gr['97.5%'], lw=1, color='k')
    ax[1].hlines(d['median'], _gr['12.5%'],
                 _gr['87.5%'], lw=2, color='k')
    ax[0].vlines(_gr['median'], _vol['2.5%'], _vol['97.5%'], lw=1, color='k')
    ax[0].vlines(_gr['median'], _vol['12.5%'], _vol['87.5%'], lw=2, color='k')
    ax[0].hlines(_vol['median'], _gr['2.5%'], _gr['97.5%'], lw=1, color='k')
    ax[0].hlines(_vol['median'], _gr['12.5%'], _gr['87.5%'], lw=2, color='k')
# plt.ylim([0, 250])
# plt.plot(d['growth_rate'], d['median'], marker='o', ms=3, color='k')

# %%
fig, ax = plt.subplots(1, 1)
for g, d in data.groupby(['carbon_source']):
    ax.hist(d['width_um'], bins=25, alpha=0.5, label=g)
ax.legend()

# %%
for g, d in data.groupby(['carbon_source']):
    plt.hist(d['width_um'], bins=20, alpha=0.5, label=g)
plt.legend()

# %%
