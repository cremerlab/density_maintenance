# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
env = data.groupby(['dataset_name', 'strain', 'condition',
                    'growth_rate_hr', 'envelope']
                   )['mass_frac'].sum().reset_index()
env = env[env['mass_frac'] <= 1]
env.sort_values(by='growth_rate_hr', inplace=True)
x = np.arange(len(env[env['envelope'] == True]['growth_rate_hr']))

fig, ax = plt.subplots(3, 1, figsize=(3, 3), sharex=True)
for i, (g, d) in enumerate(env.groupby(['dataset_name', 'condition', 'growth_rate_hr'])):
    d['mass_frac'] = d['mass_frac'] / d['mass_frac'].sum()
    _env = d[d['envelope'] == True]
    _cyt = d[d['envelope'] == False]
    ax[1].bar(x[i], _cyt['mass_frac'], color=cor['light_black'],
              edgecolor=cor['primary_black'], width=1, bottom=0)
    ax[1].bar(x[i], _env['mass_frac'],
              color=cor['light_purple'], bottom=_cyt['mass_frac'],
              width=1, edgecolor=cor['primary_black'])

data.sort_values(by='growth_rate_hr', inplace=True)
for i, (g, d) in enumerate(data[data['envelope'] == True].groupby(['growth_rate_hr', 'dataset_name', 'condition'])):
    d = d.copy()
    d['category'] = 'other'
    d['color'] = '#DACEF4'  # cor['pale_purple']
    d.loc[d['cog_letter'].isin(
        ['G', 'E', 'F', 'H', 'P', 'I']), 'category'] = 'transport & metabolism'
    d.loc[d['cog_letter'].isin(
        ['G', 'E', 'F', 'H', 'P', 'I']), 'color'] = cor['purple']
    d.loc[d['cog_letter'].isin(['N', 'T']),
          'category'] = 'motility & signaling'
    d.loc[d['cog_letter'].isin(['N', 'T']), 'color'] = cor['light_purple']
    d.loc[d['cog_letter'].isin(['C']), 'category'] = 'energy generation'
    d.loc[d['cog_letter'].isin(['C']), 'color'] = cor['primary_purple']
    d.loc[d['cog_letter'].isin(
        ['M']), 'category'] = 'cell wall/membrane biogenesis'
    d.loc[d['cog_letter'].isin(['M']), 'color'] = cor['dark_purple']
    d = d.groupby(['category', 'color'])['mass_frac'].sum().reset_index()
    d['mass_frac'] = d['mass_frac'] / d['mass_frac'].sum()
    bottom = 0
    d = d.sort_values('mass_frac', ascending=True)
    for k in ['cell wall/membrane biogenesis', 'transport & metabolism', 'energy generation', 'motility & signaling', 'other']:
        _d = d[d['category'] == k]
        ax[0].bar(x[i], _d['mass_frac'], color=_d['color'],
                  edgecolor=cor['primary_black'], width=1, bottom=bottom)
        bottom += _d['mass_frac'].values[0]


for i, (g, d) in enumerate(data[data['envelope'] == False].groupby(['growth_rate_hr', 'dataset_name', 'condition'])):
    d = d.copy()
    d['category'] = 'other'
    d['color'] = cor['light_black']
    d.loc[d['cog_letter'].isin(
        ['G', 'E', 'F', 'H', 'P', 'I']), 'category'] = 'transport & metabolism'
    d.loc[d['cog_letter'].isin(
        ['G', 'E', 'F', 'H', 'P', 'I']), 'color'] = cor['black']
    d.loc[d['cog_letter'].isin(['L']), 'category'] = 'replication'
    d.loc[d['cog_letter'].isin(['L']), 'color'] = cor['pale_black']
    d.loc[d['cog_letter'].isin(['K']), 'category'] = 'transcription'
    d.loc[d['cog_letter'].isin(['K']), 'color'] = cor['primary_black']
    d.loc[d['cog_letter'].isin(['J']), 'category'] = 'translation'
    d.loc[d['cog_letter'].isin(['J']), 'color'] = cor['dark_black']
    d = d.groupby(['category', 'color'])['mass_frac'].sum().reset_index()
    d['mass_frac'] = d['mass_frac'] / d['mass_frac'].sum()
    bottom = 0
    d = d.sort_values('mass_frac', ascending=True)
    for k in ['translation', 'transport & metabolism', 'transcription', 'other', 'replication']:
        _d = d[d['category'] == k]
        ax[2].bar(x[i], _d['mass_frac'], color=_d['color'],
                  edgecolor=cor['primary_black'], width=1, bottom=bottom)
        bottom += _d['mass_frac'].values[0]

ticks = [0, 9,  19, 29, 39, 49, 59, len(x)-1]
lam = env[env['envelope'] == True]['growth_rate_hr'].values
labels = [f'{l:0.3f}' for l in lam]
for a in ax:
    a.set_ylim([0, 1])
    a.set_yticks([0, 1])
    a.set_yticklabels([0, 1])
    a.set_xlim([-0.5, len(x)])
    a.set_xticks(x)
    a.set_xticklabels(labels, rotation=90, fontsize=3)
plt.subplots_adjust(hspace=0.1)
plt.savefig('../../figures/Fig1_bars.pdf', bbox_inches='tight')
