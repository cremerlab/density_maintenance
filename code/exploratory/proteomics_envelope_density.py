#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.empirical
import size.viz
colors, palette = size.viz.matplotlib_style()

# Load the data set
data = pd.read_csv('../../data/compiled_absolute_measurements.csv')
data.head()
periplasm = data[data['go_terms'].str.contains('GO:0042597')]
valg = periplasm[periplasm['dataset_name'] == 'Valgepea et al. 2013']
#%%
# hib = pd.DataFrame({})
# for g, d in data.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
#     tot_mass = d['fg_per_cell'].sum()
#     rmf = d[d['gene_name'].isin(['rmf'])]['fg_per_cell'].sum()
#     hib = hib.append({'growth_rate_hr': g[-1], 'dataset_name':g[0],
#                       'frac': rmf / tot_mass}, ignore_index=True)

#%%
# fig, ax = plt.subplots()
# for g, d in hib.groupby(['dataset_name']):
#     ax.plot(d['growth_rate_hr'], d['frac'] * 100, 'o', label=g)
# ax.legend()
# ax.set_xlabel('growth rate [inv hr]')
# ax.set_ylabel('mass fraction (rmf)')

#%%
subfracs = pd.DataFrame([])
for g, d in data.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    metabolic = d[d['cog_class']=='metabolism'] 
    periplasmic = metabolic[metabolic['go_terms'].str.contains('GO:0042597')]
    frac = periplasmic['fg_per_cell'].sum() / metabolic['fg_per_cell'].sum()
    subfracs = subfracs.append({'dataset_name':g[0],
                                'growth_rate_hr':g[1],
                                'condition': g[2],
                                'subfrac': frac},
                                ignore_index=True)


#%%
fig, ax = plt.subplots()
for g, d in subfracs.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['subfrac'], 'o')


#%%
mass_fracs = pd.DataFrame([])
for g, d in data.groupby(['dataset', 'dataset_name', 'growth_rate_hr', 'condition']):
    periplasm_mass = d[d['go_terms'].str.contains('GO:0042597')]['fg_per_cell'].sum()
    tot_mass = d['fg_per_cell'].sum()
    mass_fracs = mass_fracs.append({'dataset_name': g[1],
                                    'growth_rate_hr':g[2],
                                    'condition':g[3],
                                    'periplasm_mass': periplasm_mass,
                                    'total_mass': tot_mass,
                                    'mass_fraction':periplasm_mass/tot_mass},
                                    ignore_index=True)
#%%
fig, ax = plt.subplots(1, 1)
for g, d in mass_fracs.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['mass_fraction'], 'o')
ax.set_xlabel('growth rate [per hr]')
ax.set_ylabel('periplasmic protein mass fraction')

#%%

prot_frac = 350 / 500
cell_density = 500 #fg / fL 
mass_fracs['total_cell_volume'] = size.empirical.lambda2size(mass_fracs['growth_rate_hr'].values)
mass_fracs['length'] = size.empirical.lambda2length(mass_fracs['growth_rate_hr'])
mass_fracs['width'] = size.empirical.lambda2width(mass_fracs['growth_rate_hr'])
mass_fracs['envelope_volume'] = size.analytical.envelope_volume(mass_fracs['length'], mass_fracs['width'], 0.025)
mass_fracs['drymass_per_cell'] = mass_fracs['total_cell_volume'] * cell_density
mass_fracs['prot_per_cell'] = prot_frac * mass_fracs['drymass_per_cell'] 
mass_fracs['periplasmic_protein_per_cell'] = mass_fracs['prot_per_cell'] * mass_fracs['mass_fraction']

#%%
fig, ax = plt.subplots(1, 1)
for g, d in mass_fracs.groupby(['dataset_name']):
    ax.plot(d['total_mass'], d['prot_per_cell'], 'o')
    ax.set_xlabel('reported protein')
    ax.set_ylabel('our estimate')
# ax.plot(np.linspace(0, 2000), np.linspace(0, 2000), 'k--')

#%%
fig, ax = plt.subplots(1, 1)
for g, d in mass_fracs.groupby(['dataset_name']):
    ax.plot(d['periplasm_mass'] / d['envelope_volume'], (d['periplasmic_protein_per_cell']) / d['envelope_volume'] , 'o')
    ax.set_xlabel('reported protein')
    ax.set_ylabel('our estimate')
# ax.plot(np.linspace(10, 500), np.linspace(10, 500), 'k--')
ax.set_xlabel('previous density estimate (fg periplasmic protein / fL')
ax.set_ylabel('new density estimate (fg periplasmic protein / fL)')


#%%
fig, ax = plt.subplots(1, 1)
for g, d in mass_fracs.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], (d['periplasmic_protein_per_cell']) / d['envelope_volume'] , 'o')
    ax.set_xlabel('reported protein')
    ax.set_ylabel('our estimate')
# ax.plot(np.linspace(10, 500), np.linspace(10, 500), 'k--')
ax.set_xlabel('growth rate [per hr]')
ax.set_ylabel('fg periplasmic protein / fL')
ax.plot([1.2, 0.9, 0.68, 0.49], np.array([67, 72.4, 80.8, 86.5]), '*', color=colors['purple'], ms=15)

#%%
# Given the growth rate, compute the length and width
length = size.empirical.lambda2length(periplasm['growth_rate_hr'].values)
width = size.empirical.lambda2width(periplasm['growth_rate_hr'].values)
delta = 0.025
env_vol = size.analytical.envelope_volume(length, width, delta)
periplasm['envelope_volume'] = env_vol

# %%
# Compute the total protein density in periplasm
density_df = pd.DataFrame({})
for g, d in periplasm.groupby(['dataset_name', 'growth_rate_hr', 'condition']):
    vol = d['envelope_volume'].values[0]
    tot_mass = d['fg_per_cell'].sum()
    tot_density = tot_mass / vol
    density_df = density_df.append({'growth_rate_hr':g[1],
                                    'dataset_name':g[0],
                                    'condition':g[-1],
                                    'tot_mass': tot_mass,
                                    'envelope_vol': vol,
                                    'density': tot_density},
                                    ignore_index=True)
                                    



density_df.to_csv('../../data/periplasmic_protein_density.csv', index=False)

# %%
fig, ax = plt.subplots(1,1)
for g, d in density_df.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['density']/1000, 'o', label=g, ms=4)

ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('periplasmic protein density [pg / fL]')
ax.legend()
ax.set_ylim([0, 0.4])
ax.set_xlim([0, 2])
plt.savefig('./periplasmic_density.pdf', bbox_inches='tight')
# %%
fig, ax = plt.subplots(1,1)
for g, d in density_df.groupby(['dataset_name']):
    ax.plot(d['growth_rate_hr'], d['tot_mass'], 'o', label=g, ms=4)

ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('periplasmic protein mass per cell [fg]')
ax.legend()
ax.set_ylim([0, 30])
ax.set_xlim([0, 2])
plt.savefig('./periplasmic_mass_per_cell.pdf', bbox_inches='tight')
# %%
