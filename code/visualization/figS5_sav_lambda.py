#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load 
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')

# Directly calculate the emprical kappa given by Eq. 8 in main text
W_PERI = 0.025
BETA = 2.19
data['empirical_kappa'] = (2 * data['surface_area_um2'] / (data['volume_fL'] - W_PERI * data['surface_area_um2'])) * (data['phi_cyto'] + BETA * data['phi_rib']) / data['phi_mem']

#%% Plot the break of the sav and growth rate
# Set figure canvas
fig, ax = plt.subplots(1, 2, figsize=(4, 1.75))

# Plot the wildtype data
wt = data[data['strain']=='wildtype']
ax[0].plot(wt['growth_rate_hr'], wt['surface_to_volume_inv_um'], 'o',
        markerfacecolor=cor['pale_black'], ms=4,
        markeredgecolor=cor['pale_black'], alpha=0.5)
ax[1].plot(wt['growth_rate_hr'], wt['empirical_kappa'], 'o',
        markerfacecolor=cor['pale_black'], ms=4,
        markeredgecolor=cor['pale_black'], alpha=0.5)


# Define colors and markers for the perturbations.
inducer_concs = {0:'pale_', 2:'primary_', 4:'', 100:''}
strain_colors = {'relA':  'red',
                 'meshI': 'gold'} 
markers = {'glucose': 'D', 'glucoseCAA': 's'}

for g, d in data[data['strain']!='wildtype'].groupby(['carbon_source', 'strain', 'inducer_conc']):
    if g[1] == 'meshI':
        zorder=1000
    else:
        zorder=10

        
    fmt = {'marker':markers[g[0]],
           'markeredgecolor':cor[f'dark_{strain_colors[g[1]]}'],
           'color':cor[f'{inducer_concs[g[2]]}{strain_colors[g[1]]}'],
           'alpha':0.85,
           'markersize':4,
           'linestyle':'none',
           'zorder':zorder}
    ax[0].plot(d['growth_rate_hr'], d['surface_to_volume_inv_um'],
            **fmt)
    ax[1].plot(d['growth_rate_hr'], d['empirical_kappa'],
            **fmt)
# Add
# Add context
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
ax[1].set_ylabel('cytoplasm-membrane\ndensity ratio [µm$^{-1}$]', fontsize=6)
ax[0].set_xlim([0.25, 1.5])
ax[1].set_xlim([0.25, 1.5])
ax[0].set_ylim([4.5, 8])
ax[1].set_ylim([50, 250])
plt.savefig('./plots/figS5_sav_lambda_plot.pdf', bbox_inches='tight')



