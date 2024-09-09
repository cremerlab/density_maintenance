#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load our data 
data = pd.read_csv('../../data/collated/collated_mass_fractions.csv')

# Define the types
cog_class = {'info': [['J', 'A', 'K', 'L','B'], 'black'],
             'signaling': [['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'U', 'O'], 'red'],
             'metabolism': [['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'], 'gold'],
             'other': [['R', 'S', 'X'], 'purple']}

data['color'] = 'black' 
data['class'] = 'other'
for k, v in cog_class.items():
    data.loc[data['cog_letter'].isin(v[0]), 'class'] = k
    data.loc[data['cog_letter'].isin(v[0]), 'color'] = v[1]

# Define the important localization
data['compartment'] = 'envelope'
data.loc[data['localization']=='cytoplasm', 'compartment'] = 'cytoplasm'
data.loc[data['replicate'].isnull(), 'replicate'] = 0

# Group and sum for plotting
cog_data = data.groupby(['source', 'condition', 'growth_rate_hr', 'compartment', 'color', 'class', 
                    'replicate', 'strain'])[['mass_frac']].sum().reset_index()

fig, ax = plt.subplots(2, 1, figsize=(3, 3))
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('composition of\ncytoplasmic proteome', fontsize=6)
ax[1].set_ylabel('composition of\nenvelope proteome', fontsize=6)

axes = {'cytoplasm': ax[0], 'envelope': ax[1]}
for g, d in cog_data.groupby(['compartment','source', 'strain', 'replicate', 'class', 'color']):
    fmt = size.viz.style_point(g[1], alpha=0.5)
    fmt['color'] = cor[f'primary_{g[-1]}']
    if g[1] == 'This Study':
        fmt['color'] = cor[f'pale_{g[-1]}'] 
        fmt['markeredgewidth'] = 0.75 
        fmt['markeredgecolor'] = cor[f'dark_{g[-1]}']
        fmt['alpha'] = 1
    axes[g[0]].plot(d['growth_rate_hr'], d['mass_frac'], **fmt)
ax[0].set_ylim([-0.05, 0.6])
ax[1].set_ylim([-0.01, 0.18])

plt.savefig('./plots/fig1_compartment_COG_trends.pdf', bbox_inches='tight')


#%%
# Map the ribosomal proteins
ribo_prots = ['rrsa', 'rpsa', 'rpsb', 'rpsc', 'rpd', 'rpse', 'rpsf', 'rpsg', 'rpsh', 'rpsi', 'rpsj', 'rpsk',
              'rpsl', 'rpsm', 'rpsn', 'rpso', 'rpsp', 'rpsq', 'rpsr', 'rpst', 'rpsu', 'rrla', 'rrfa', 'rpla', 'rplb', 'rplc', 'rpld', 'rple', 'rplf', 'rplj', 'rpll', 'rpli',
              'rplk', 'rplm', 'rpln', 'rplo', 'rplop', 'rplq', 'rplr', 'rpls', 'rplt', 'rplu', 'rplv', 'rplw', 'rplx', 'rply', 'rpma',
              'rpmb', 'rpmc', 'rpmd', 'rpme', 'rpmf', 'rmpg', 'rpmh', 'rpmi', 'rpmj']

localizations = data.copy()
mems = ['inner membrane', 'outer membrane', 'membrane related']
localizations.loc[localizations['localization']=='cytoplasm', 'localization'] = 'cytoplasmic'
localizations.loc[localizations['localization']=='periplasm', 'localization'] = 'periplasmic'
localizations.loc[localizations['localization'].isin(mems), 'localization'] = 'membrane'
localizations = localizations.groupby(['source', 'strain', 'condition', 'growth_rate_hr', 'replicate', 'localization'])[['mass_frac']].sum().reset_index()

fig, ax = plt.subplots(1, 3, figsize=(6, 1.5))
# for g, d in data.groupby(['source', 'strain', 'condition', 'growth_rate_hr', 'replicate']):
#     ribo_frac = d[d['gene_name'].str.lower().isin(ribo_prots)]['mass_frac'].sum()
#     fmt = size.viz.style_point(g[0], alpha=0.5)
#     if g[0] == 'This Study':
#         zorder = 1000
#     else:
#         zorder = 10
#     ax[0].plot(g[-2], ribo_frac, **fmt, zorder=zorder)
for g, d in localizations.groupby(['source']):
    cyto = d[d['localization']=='cytoplasmic']
    peri = d[d['localization']=='periplasmic']
    memb = d[d['localization']=='membrane']
    fmt = size.viz.style_point(g[0], alpha=0.5)
    ax[0].plot(cyto['growth_rate_hr'], cyto['mass_frac'], **fmt)
    if g[0] == 'This Study':
        fmt['color'] = cor['pale_purple']
        fmt['markeredgecolor'] = cor['dark_purple']
    else:
        fmt['color'] = cor['primary_purple']
    ax[1].plot(peri['growth_rate_hr'], peri['mass_frac'], **fmt)

    if g[0] == 'This Study':
        fmt['color'] = cor['pale_blue']
        fmt['markeredgecolor'] = cor['dark_blue']
    else:
        fmt['color'] = cor['primary_blue']
 
    ax[2].plot(memb['growth_rate_hr'], memb['mass_frac'], **fmt)

ax[0].set_ylim([0.5, 1.0])
ax[1].set_ylim([0, 0.15])
ax[2].set_ylim([0, 0.2])
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$\phi_{cyt}$\ncytoplasmic\nproteome fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)
ax[2].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/fig1_compartment_localization.pdf', bbox_inches='tight')
#%%