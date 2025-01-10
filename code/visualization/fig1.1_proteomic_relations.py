#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
import scipy.stats
cor, pal = size.viz.matplotlib_style()

# Load our data 
data = pd.read_csv('../../data/collated/merged_mass_spectrometry.csv')

# Define the types
cog_class = {'info': [['J', 'A', 'K', 'L','B'], 'blue'],
             'signaling': [['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'U', 'O'], 'green'],
             'metabolism': [['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'], 'purple'],
             'other': [['R', 'S', 'X'], 'black']}

data['color'] = 'black' 
data['class'] = 'other'
for k, v in cog_class.items():
    data.loc[data['cog_letter'].isin(v[0]), 'class'] = k
    data.loc[data['cog_letter'].isin(v[0]), 'color'] = v[1]

#%%
# # Define the important localization
data['compartment'] = 'envelope'
data.loc[data['localization']=='phi_cyto', 'compartment'] = 'cytoplasm'
data.loc[data['replicate'].isnull(), 'replicate'] = 0

# Group and sum for plotting
cog_data = data.groupby(['source', 'carbon_source', 'strain', 'growth_rate_hr', 'compartment', 'color', 'class', 
                    'replicate',])[['mass_frac']].sum().reset_index()

fig, ax = plt.subplots(2, 1, figsize=(3, 3))
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('composition of\ncytoplasmic proteome', fontsize=6)
ax[1].set_ylabel('composition of\nenvelope proteome', fontsize=6)

axes = {'cytoplasm': ax[0], 'envelope': ax[1]}
for g, d in cog_data.groupby(['compartment' ,'source', 'strain', 'replicate', 'class', 'color']):
    fmt = size.viz.style_point(g[1], alpha=0.5)
    fmt['color'] = cor[f'primary_{g[-1]}']
    if g[1] == 'This Study':
        fmt['color'] = cor[f'pale_{g[-1]}'] 
        fmt['markeredgewidth'] = 0.75 
        fmt['markeredgecolor'] = cor[f'dark_{g[-1]}']
        fmt['alpha'] = 1

    # if (g[1] == 'This Study') & (g[0] == 'envelope'):
        # break
    axes[g[0]].plot(d['growth_rate_hr'], d['mass_frac'], **fmt)
# ax[0].set_ylim([-0.05, 0.6])
# ax[1].set_ylim([-0.01, 0.18])

# plt.savefig('./plots/fig1_compartment_COG_trends.pdf', bbox_inches='tight')


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

#%%
# Plot the 'envelope growth laws'
fig, ax = plt.subplots(1, 3, figsize=(6, 1.5))
for g, d in data.groupby(['source', 'strain', 'condition', 'growth_rate_hr', 'replicate']):
    ribo_frac = d[d['gene_name'].str.lower().isin(ribo_prots)]['mass_frac'].sum()
    fmt = size.viz.style_point(g[0], alpha=0.4)
    if g[0] == 'This Study':
        zorder = 1000
        fmt['color'] = cor['pale_gold']
        fmt['markeredgecolor'] = cor['dark_gold']
    else:
        fmt['color'] = cor['primary_gold']
        fmt['markeredgecolor'] = cor['primary_black']
        zorder = 10
    ax[0].plot(g[-2], ribo_frac, **fmt, zorder=zorder)


for g, d in localizations.groupby(['source']):
    # cyto = d[d['localization']=='cytoplasmic']
    peri = d[d['localization']=='periplasmic']
    memb = d[d['localization']=='membrane']
    fmt = size.viz.style_point(g[0], alpha=0.4)
    # ax[0].plot(cyto['growth_rate_hr'], cyto['mass_frac'], **fmt)
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

# Set the growth rate range for empirical fits
lam_range = np.linspace(0, 2.5, 200)

# Empirical fits to ribosomal data
phi_rib_data = data[(data['source']=='This Study') & (data['gene_name'].str.lower().isin(ribo_prots))]
phi_rib_data = phi_rib_data.groupby(['strain', 'condition', 'growth_rate_hr', 'replicate'])['mass_frac'].sum().reset_index() 
phi_rib_popt = scipy.stats.linregress(phi_rib_data['growth_rate_hr'], phi_rib_data['mass_frac'])
phi_rib_fit = phi_rib_popt[1] + phi_rib_popt[0] * lam_range
ax[0].plot(lam_range, phi_rib_fit, '--', lw=1,color=cor['dark_gold'], zorder=1000)

# Empirical fits to periplasmic data
phi_peri_data = localizations[(localizations['source']=='This Study') & (localizations['localization']=='periplasmic')]
phi_peri_popt = scipy.stats.linregress(phi_peri_data['growth_rate_hr'], np.log(phi_peri_data['mass_frac']))
phi_peri_fit = np.exp(phi_peri_popt[1] + phi_peri_popt[0] * lam_range)
ax[1].plot(lam_range, phi_peri_fit, '--', lw=1,color=cor['dark_purple'], zorder=1000)

# Empirical fit to membrane data
phi_mem_data = localizations[(localizations['source']=='This Study') & (localizations['localization']=='membrane')]
phi_mem_popt = scipy.stats.linregress(phi_mem_data['growth_rate_hr'], phi_mem_data['mass_frac'])
phi_mem_fit = phi_mem_popt[1] + phi_mem_popt[0] * lam_range
ax[2].plot(lam_range, phi_mem_fit, '--', lw=1,color=cor['dark_blue'], zorder=1000)

# ax[0].set_ylim([0.5, 1.0])
ax[1].set_ylim([0, 0.15])
ax[2].set_ylim([0, 0.2])
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$\phi_{cyt}$\ncytoplasmic\nproteome fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)
ax[2].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/fig1_composition_empirics.pdf', bbox_inches='tight')

#%%
# Generate an easy legend
fig, ax = plt.subplots(1,1)

for g, d in data.groupby(['source']):
    fmt = size.viz.style_point(g[0])
    ax.plot([], [], **fmt)
ax.legend()
plt.savefig('./plots/fig1_legend.pdf', bbox_inches='tight')
