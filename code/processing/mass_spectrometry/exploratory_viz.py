#%%
import scipy.stats
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 
import size.viz
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()

# Load the literature proteomics data
lit_data = pd.read_csv('../../../data/literature/collated_mass_fractions_empirics.csv')

# Load our data set and compute the theory.
data = pd.read_csv('./total_collated_data.csv')
rel_data = pd.read_csv('./total_collated_data_relative.csv')
def theory(phi_mem, phi_peri, phi_rib, kappa=110, beta=0.4558):
    numer = phi_mem * kappa
    denom = 2 * (1 + phi_rib/beta - phi_mem - phi_peri)
    return numer / denom

# Load the protein per cell data to calculate the densities
prot_per_cell = pd.read_csv('../../../data/literature/collated_protein_per_cell.csv')
drymass_density = 287
popt = scipy.stats.linregress(prot_per_cell['growth_rate_hr'], np.log(prot_per_cell['fg_protein_per_cell']))
data['prot_per_cell'] = np.exp(popt[0] * data['growth_rate_hr'] + popt[1]) 
data['membrane_mass'] = data['phi_mem'] * data['prot_per_cell']
data['membrane_density'] = data['membrane_mass'] / (2 * data['surface_area'])
data['kappa'] = drymass_density / data['membrane_density']
data['pred_SAV'] = theory(data['phi_mem'], data['phi_peri'], data['phi_rib'], data['kappa'])

#%%
fig, ax = plt.subplots(1, 3, figsize=(6,2))
wt_data = data[data['strain']=='wildtype']
wt_rel_data = rel_data[rel_data['strain']=='wildtype']


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

relA= data[data['strain']=='relA']
relA.loc[relA['inducer_conc']==0, 'color'] = cor['pale_red'] 
relA.loc[relA['inducer_conc']==2, 'color'] = cor['primary_red'] 
relA.loc[relA['inducer_conc']==4, 'color'] = cor['dark_red'] 

for g, d in relA.groupby('color'):
    style['markerfacecolor'] = g
    ax[0].plot(d['growth_rate_hr'], d['phi_rib'], 'o',**style)
    ax[1].plot(d['growth_rate_hr'], d['phi_mem'], 'o',**style)
    ax[2].plot(d['growth_rate_hr'], d['phi_peri'], 'o',**style)


meshI = data[data['strain']=='meshI']
meshI.loc[meshI['inducer_conc']==0, 'color'] = cor['pale_green'] 
meshI.loc[meshI['inducer_conc']==100, 'color'] = cor['dark_green'] 


for g, d in meshI.groupby('color'):
    style['markerfacecolor'] = g
    ax[0].plot(d['growth_rate_hr'], d['phi_rib'], 'o',**style)
    ax[1].plot(d['growth_rate_hr'], d['phi_mem'], 'o',**style)
    ax[2].plot(d['growth_rate_hr'], d['phi_peri'], 'o',**style)



lacZ = data[data['strain']=='lacZ']
lacZ.loc[lacZ['inducer_conc']==0, 'color'] = cor['pale_blue'] 
lacZ.loc[lacZ['inducer_conc']==1, 'color'] = cor['light_blue'] 
lacZ.loc[lacZ['inducer_conc']==3, 'color'] = cor['primary_blue'] 
lacZ.loc[lacZ['inducer_conc']==5, 'color'] = cor['blue'] 

for g, d in lacZ.groupby('color'):
    style['markerfacecolor'] = g
    ax[0].plot(d['growth_rate_hr'], d['phi_rib'], 'o',**style)
    ax[1].plot(d['growth_rate_hr'], d['phi_mem'], 'o',**style)
    ax[2].plot(d['growth_rate_hr'], d['phi_peri'], 'o',**style)




for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$\phi_{rib}$\n ribosomal mass fraction', fontsize=6)
ax[1].set_ylabel('$\phi_{mem}$\n membrane mass fraction', fontsize=6)
ax[2].set_ylabel('$\phi_{peri}$\n periplasmic mass fraction', fontsize=6)
ax[0].set_ylim([0, 0.3])
ax[1].set_ylim([0.05, 0.2])
ax[2].set_ylim([0, 0.15])
plt.tight_layout()
# ax[0].plot(lit_phi_rib[''])
# ax[0].legend()

#%%
lit_vol_data = pd.read_csv('../../../data/literature/collated_literature_volume_data.csv')
lit_size_data = pd.read_csv('../../../data/literature/collated_literature_size_data.csv')
lit_size_data['alpha'] = lit_size_data['length_um'] / lit_size_data['width_um']
fig, ax = plt.subplots(2,2, figsize=(4, 3))
ax = ax.ravel()
for g, d in lit_vol_data.groupby('source'):
    style = size.viz.style_point(g) 
    ax[-1].plot(d['growth_rate_hr'], d['volume_um3'], **style)
    
for g, d in lit_size_data.groupby('source'):
    style = size.viz.style_point(g) 
    ax[0].plot(d['growth_rate_hr'], d['width_um'], **style)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], **style)
    ax[2].plot(d['growth_rate_hr'], d['alpha'], **style)

style = {"ms": 4, "marker": 'o', 'markerfacecolor': 'w', 'markeredgecolor': cor['primary_black'],
         "markeredgewidth":1, 'linestyle':'none'}
ax[0].plot(wt_data['growth_rate_hr'], wt_data['width_median'], **style)
ax[1].plot(wt_data['growth_rate_hr'], wt_data['length'], **style)
ax[2].plot(wt_data['growth_rate_hr'], wt_data['length'] / wt_data['width_median'], **style)
ax[3].plot(wt_data['growth_rate_hr'], wt_data['volume'], **style)

ax[0].set_ylim([0.4, 1.2])
ax[1].set_ylim([1, 4])
ax[2].set_ylim([1, 5])
ax[3].set_ylim([0, 5])

for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('width [$\mu$m]', fontsize=6)
ax[1].set_ylabel('length [$\mu$m]', fontsize=6)
ax[2].set_ylabel('aspect ratio', fontsize=6)
ax[3].set_ylabel('volume [fL]', fontsize=6)
plt.tight_layout()


#%%
plt.plot([4, 10], [4,10], 'k:')
plt.plot(wt_data['pred_SAV'], wt_data['surface_to_volume'], 'o', alpha=0.5, label='wildtype')
relA.sort_values('inducer_conc', inplace=True)
plt.plot(relA['pred_SAV'], relA['surface_to_volume'], '-', lw=0.75, color=cor['red'])
for g, d in relA.groupby('inducer_conc', sort='ascending'):
    plt.plot(d['pred_SAV'], d['surface_to_volume'], 'o', color=d.color.values[0], 
    label=f'relA {g} ng/mL DOX') 
meshI.sort_values('inducer_conc', inplace=True)
plt.plot(meshI['pred_SAV'], meshI['surface_to_volume'], '-', lw=0.75, color=cor['green'])
for g, d in meshI.groupby('inducer_conc', sort='ascending'):
    plt.plot(d['pred_SAV'], d['surface_to_volume'], '-o', color=d.color.values[0],
    label=f'meshI {g} mM IPTG')  

lacZ.sort_values('inducer_conc', inplace=True)
plt.plot(lacZ['pred_SAV'], lacZ['surface_to_volume'], '-', lw=0.75, color=cor['blue'])
for g, d in lacZ.groupby('inducer_conc', sort='ascending'):
    plt.plot(d['pred_SAV'], d['surface_to_volume'], '-o', color=d.color.values[0],
    label=f'lacZ {g} ng/mL CTC') 

plt.legend()
plt.xlabel('predicted SA/V [µm$^{-1}$]', fontsize=6)
plt.ylabel('measured SA/V [µm$^{-1}$]', fontsize=6)

#%%
# Load the empirical density data
emp_density = pd.read_csv('../../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
kappa = emp_density[emp_density['quantity']=='ms_rho_mem']

fig, ax = plt.subplots(1,1, figsize=(3, 2))
ax.set_ylim([0, 6])
for g, d in kappa.groupby('source'):
    style = size.viz.style_point(g) 
    plt.plot(d['growth_rate_hr'], d['median_value'], **style)

ax.plot(wt_data['growth_rate_hr'], wt_data['membrane_density'],'o',
        markerfacecolor='w', markeredgecolor=cor['primary_black'], markeredgewidth=1,
        label='wildtype')

for g, d in relA.groupby('color'):
    ax.plot(d['growth_rate_hr'], d['membrane_density'],'o',
              markerfacecolor=g, markeredgewidth=1,
             markeredgecolor='k', label='relA')

for g, d in meshI.groupby('color'):
    ax.plot(d['growth_rate_hr'], d['membrane_density'],'o',
              markerfacecolor=g, markeredgewidth=1,
             markeredgecolor='k', label='meshI')

for g, d in lacZ.groupby('color'):
    ax.plot(d['growth_rate_hr'], d['membrane_density'],'o',
              markerfacecolor=g, markeredgewidth=1,
             markeredgecolor='k', label='lacZ')
