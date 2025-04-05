#%%
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the mass spec data and group by cog letter/compartment
ms_data = pd.read_csv('../../data/collated/experimental_mass_spectrometry.csv')
ms_data = ms_data[ms_data['strain']=='wildtype']

# Apply higher-order grouping
rev_mapper = {'metabolism': ['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'],
          'signaling': ['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'],
          'info': ['J', 'A', 'K', 'L', 'B'],
          'other': ['R', 'S', 'X']}
mapper = {_v:k for k, v in rev_mapper.items() for _v in v}
ms_data['cog_group'] = [mapper[ell] for ell in ms_data['cog_letter'].values]

# Group and aggregate by cog group
ms_data = ms_data.groupby(['date', 'carbon_source', 'replicate', 'growth_rate_hr', 
                           'localization', 'cog_group'])['mass_frac'].sum().reset_index()


#%% Instantiate the figure canvas
fig, ax = plt.subplots(1, 3, figsize=(6, 2), sharex=True)
axes = {'psi_cyto': ax[0], 'psi_peri': ax[1], 'psi_mem':ax[2]}

# Set colors for different compartments
pal = sns.color_palette('husl', 4)
mapper = {'info': pal[0], 'metabolism': pal[1], 'other': pal[2], 'signaling':pal[3]}

# Style the marker
fmt = size.viz.style_point('This Study')
for g, d in ms_data.groupby(['localization', 'date', 'carbon_source', 'replicate']):
    if g[0] == 'psi_ext':
        continue
    # Compute the relative mass fraction change in the compartment
    tot_mass = d['mass_frac'].sum()
    d['rel_mass_frac'] = d['mass_frac'] / tot_mass

    # Plot
    for _g, _d in d.groupby('cog_group'):
        fmt['markerfacecolor'] = mapper[_g]
        fmt['markeredgewidth'] = 0.75
        axes[g[0]].plot(_d['growth_rate_hr'], _d['rel_mass_frac'], **fmt)

# Add context.
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
    a.set_ylim([-0.05, 1])

ax[0].set_ylabel('cytoplasm proteome\nfunctional composition', fontsize=6)
ax[1].set_ylabel('periplasm proteome\nfunctional composition', fontsize=6)
ax[2].set_ylabel('membrane proteome\nfunctional composition', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/figS3_cog_compartments_plots.pdf', bbox_inches='tight')