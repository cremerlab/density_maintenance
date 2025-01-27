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
cog_class = {'info': [['J', 'A', 'K', 'L','B'], 'goldenrod'],
             'signaling': [['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'U', 'O'], 'yellowgreen'],
             'metabolism': [['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'], 'firebrick'],
             'other': [['R', 'S', 'X'], 'gray']}

data['color'] = 'gray' 
data['class'] = 'other'
for k, v in cog_class.items():
    data.loc[data['cog_letter'].isin(v[0]), 'class'] = k
    data.loc[data['cog_letter'].isin(v[0]), 'color'] = v[1]

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
for g, d in cog_data.groupby(['compartment', 'source', 'strain', 'carbon_source', 'replicate']):
    # Compute the total mass of the compartment
    compartment_mass = d['mass_frac'].sum()

    for _g, _d in d.groupby(['class', 'color']):
        fmt = size.viz.style_point(g[1], alpha=0.4)
        fmt['color'] = _g[-1]
        if g[1] == 'This Study':
            # fmt['color'] = cor[f'pale_{_g[-1]}'] 
            fmt['markeredgewidth'] = 0.75 
            # fmt['markeredgecolor'] = cor[f'dark_{_g[-1]}']
            fmt['alpha'] = 0.75 
        axes[g[0]].plot(_d['growth_rate_hr'], _d['mass_frac'] / compartment_mass, **fmt)
plt.savefig('./plots/fig1_compartment_COG_trends.pdf', bbox_inches='tight')

#%%
# Compute the total proteome mass fractions in each compartment
localization = data.groupby(['source', 'localization', 'strain', 'carbon_source', 
                             'growth_rate_hr', 'replicate'])['mass_frac'].sum().reset_index()

# Define the figure canvas and set the axes / colors
fig, ax = plt.subplots(1, 3, figsize=(6, 1.5), sharex=True)
axes = {'phi_cyto': [ax[0], 'black'], 
        'phi_peri': [ax[1], 'purple'], 
        'phi_mem': [ax[2], 'blue']}

# Define the growth rate range over which to compute the simplistic fit
lam_range = np.linspace(0, 2.5, 100)
res = {}  # Result dictionary for the fit. 

# Iterate through each source and compartment and plot
for g, d in localization.groupby(['source', 'localization']):
    fmt = size.viz.style_point(g[0], alpha=0.2)
    fmt['color'] = cor[f'primary_{axes[g[1]][1]}']

    # Set specific styling if it's this study
    if g[0] == 'This Study':
        fmt['markeredgecolor'] = cor[f'primary_{axes[g[1]][1]}']
        fmt['markerfacecolor'] = 'white'
        fmt['markeredgewidth'] = 1
        fmt['alpha'] = 1
        
     
    _ax = axes[g[1]][0]
    _ax.plot(d['growth_rate_hr'], d['mass_frac'], **fmt)


    # If source is this study, perform the simplistic fits and plot
    if g[0] == 'This Study':
        if g[1] == 'phi_peri':
            popt = scipy.stats.linregress(d['growth_rate_hr'], np.log(d['mass_frac']))
            fit = np.exp(popt[1] + popt[0] * lam_range)
        else:
            popt = scipy.stats.linregress(d['growth_rate_hr'], d['mass_frac'])
            fit = popt[1] + popt[0] * lam_range 
        res[g[1]] = popt
        _ax.plot(lam_range, fit, '--', lw=1.5, color=fmt['color'])



# Set context
ax[0].set_ylabel('$\phi_{cyt}$\ncytoplasmic\nproteome fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)
ax[2].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0].set_ylim([0.6, 0.9])
ax[1].set_ylim([0, 0.15])
ax[1].set_yticks([0, 0.05, 0.1, 0.15])
ax[2].set_yticks([0, 0.05, 0.1, 0.15, 0.2])
ax[2].set_ylim([0, 0.2])
plt.savefig('./plots/fig1_localization_trends.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1,1, figsize=(1.5, 1.5))
for g, d in localization.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    fmt['color'] = cor[f'primary_black']

    # Set specific styling if it's this study
    if g == 'This Study':
        fmt['markeredgecolor'] = cor[f'dark_black']
        fmt['markerfacecolor'] = 'white'
        fmt['markeredgewidth'] = 1
        fmt['alpha'] = 1 

    peri = d[d['localization']=='phi_peri']
    cyto = d[d['localization']=='phi_cyto']
    ax.plot(peri['growth_rate_hr'], peri['mass_frac'].values + cyto['mass_frac'].values,
            **fmt)
    ax.set_ylim([0.6, 0.95])

    # If our data, compute the linear fit. 
    if g == 'This Study':
        popt = scipy.stats.linregress(peri['growth_rate_hr'], peri['mass_frac'].values + cyto['mass_frac'].values)
        fit = popt[1] + popt[0] * lam_range
        ax.plot(lam_range, fit, 'k--', lw=1.5)


ax.legend(bbox_to_anchor=(1,1))
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('$\phi_{cyto} + \phi_{peri}$\ncytoplasmic + periplasmic \nproteome allocation', fontsize=6)
plt.savefig('./plots/figS1_cytoplasmic_periplasmic_sum.pdf')