#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('../../data/compiled_measurements.csv')

kappa = 134
data['pred_SAV'] = kappa * data['phi_mem'] / (2 * (1 + (data['phi_rib']/0.4558) - data['phi_mem'] - data['phi_peri']))

fig, ax = plt.subplots(1,1, figsize=(4,4))
wt = data[data['strain']=='wildtype']

ax.plot(wt['sav_inv_um'], wt['pred_SAV'], 'o',
        markeredgewidth=1, markeredgecolor='k',
        markerfacecolor='w', ms=4)

carbon_sources = {'glucose': 'o',
                  'glucoseCAA': 'D'}
ind_colors = {100: 'dark_', 0:'pale_', 2:'', 4:'dark_'} 
strains = {'relA': 'red', 'meshI': 'green'}
for g, d in data[data['strain']!='wildtype'].groupby(['carbon_source', 'strain', 'inducer_conc']):
    ax.plot(d['sav_inv_um'], d['pred_SAV'],
            marker=carbon_sources[g[0]],
            markeredgecolor=cor['primary_black'],
            markeredgewidth=1, 
            markersize=5,
            markerfacecolor=cor[f'{ind_colors[g[2]]}{strains[g[1]]}'],
            label=f'carbon: {g[0]}\nstrain: {g[1]}\nind. level:{g[1]}',
            linestyle='none')

ax.plot([4, 8], [4, 8], 'k--')

ax.set_xlabel('measured surface-to-volume [µm$^{-1}$]', fontsize=6)
ax.set_ylabel('predicted surface-to-volume [µm$^{-1}$]', fontsize=6)

#%%
fig, ax = plt.subplots(1,1)
for g, d in data[data['strain']!='wildtype'].groupby(['carbon_source', 'strain', 'inducer_conc']):
    ax.plot(d['growth_rate_hr'], d['phi_mem'],
            marker=carbon_sources[g[0]],
            markeredgecolor=cor['primary_black'],
            markeredgewidth=1, 
            markersize=5,
            markerfacecolor=cor[f'{ind_colors[g[2]]}{strains[g[1]]}'],
            label=f'carbon: {g[0]}\nstrain: {g[1]}\nind. level:{g[1]}',
            linestyle='none')
ax.plot(wt['growth_rate_hr'], wt['phi_mem'], 'o',
        markeredgewidth=1, markeredgecolor='k',
        markerfacecolor='w', ms=4)

ax.set_ylabel('$\phi_{mem}$')
ax.set_xlabel('growth rate [hr$^{-1}$]')

ax.set_ylim([0, 0.2])

#%%
fig, ax = plt.subplots(1,1)
for g, d in data[data['strain']!='wildtype'].groupby(['carbon_source', 'strain', 'inducer_conc']):
    ax.plot(d['growth_rate_hr'], d['phi_rib'],
            marker=carbon_sources[g[0]],
            markeredgecolor=cor['primary_black'],
            markeredgewidth=1, 
            markersize=5,
            markerfacecolor=cor[f'{ind_colors[g[2]]}{strains[g[1]]}'],
            label=f'carbon: {g[0]}\nstrain: {g[1]}\nind. level:{g[1]}',
            linestyle='none')
ax.plot(wt['growth_rate_hr'], wt['phi_rib'], 'o',
        markeredgewidth=1, markeredgecolor='k',
        markerfacecolor='w', ms=4)

ax.set_ylabel('$\phi_{rib}$')
ax.set_xlabel('growth rate [hr$^{-1}$]')
# ax.set_xlim([0, 1.5])
# ax.set_ylim([0, 0.2])