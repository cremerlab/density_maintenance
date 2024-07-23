#%%
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 
import size.viz
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

# Load the literature proteomics data
lit_data = pd.read_csv('../../../data/literature/collated_mass_fractions_empirics.csv')

# Load this dataset
data = pd.read_csv('./total_collated_data.csv')
wt_data = data[data['strain']=='wildtype']

# lit_phi_rib = pd.read_csv('../../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
#%%

fig, ax = plt.subplots(1, 3, figsize=(6,2))

for g, d in lit_data.groupby('dataset_name'):
    lit_phi_rib = d[d['localization'] == 'ribosomal sector']
    lit_phi_mem = d[d['localization'] == 'membrane']
    lit_phi_peri = d[d['localization'] == 'periplasm']
    pt = mapper[g]
    style = {'ms': 4, 'marker': pt['m'], 'markerfacecolor': pt['c'], 'markeredgecolor': 'k', 'label': g,
             'alpha':0.5}
    for i, _d in enumerate([lit_phi_rib, lit_phi_mem, lit_phi_peri]):
        ax[i].plot(_d['growth_rate_hr'], _d['mass_frac'], linestyle='none', **style)

style = {"ms": 4, "marker": 'o', 'markerfacecolor': 'w', 'markeredgecolor': cor['primary_black'],
         "markeredgewidth":1}
ax[0].plot(wt_data['growth_rate_hr'], wt_data['phi_rib'], 'o',**style)
ax[1].plot(wt_data['growth_rate_hr'], wt_data['phi_mem'], 'o',**style)
ax[2].plot(wt_data['growth_rate_hr'], wt_data['phi_peri'], 'o',**style)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylim([0, 0.3])
ax[1].set_ylim([0, 0.2])
ax[2].set_ylim([0, 0.15])

# ax[0].plot(lit_phi_rib[''])
# ax[0].legend()
