#%%

import numpy as np 
import pandas as pd 
import size.viz 
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
cor, pal = size.viz.matplotlib_style()

# Load our data
data = pd.read_csv('../../processing/mass_spectrometry/total_collated_data.csv')

#Load collated ms data
ms_data = pd.read_csv('../../../data/literature/collated_mass_fractions_empirics.csv')

#Load the literature data
prot_per_cell = pd.read_csv('../../../data/literature/collated_protein_per_cell.csv')
dna_data = pd.read_csv('../../../data/literature/collated_dna_protein_ratio.csv')
size_data = pd.read_csv('../../../data/literature/collated_literature_size_data.csv')
vol_data = pd.read_csv('../../../data/literature/collated_literature_volume_data.csv')

#%%

# Load the cmdstan model
model = cmdstanpy.CmdStanModel(stan_file='./revised_density_ratio.stan')

#%%
lam_range = np.linspace(0.1, 2.5, 200)
# Assemble the data dictionary
data_dict = {
    'N_prot':len(prot_per_cell),
    'N_SA': len(size_data),
    'N_vol': len(vol_data),
    'N_DNA': len(dna_data),
    'N_lit_ms': len(ms_data[ms_data['localization']=='envelope']),
    'N_ms': len(data),
    'N_scan': len(lam_range),
    'scan_lam': lam_range,

    'protein_per_cell': prot_per_cell['fg_protein_per_cell'].values,
    'protein_per_cell_lambda': prot_per_cell['growth_rate_hr'].values,
    'SA': size_data['surface_area_um2'].values,
    'SA_lambda': size_data['growth_rate_hr'].values,
    'cell_volume': vol_data['volume_um3'].values,
    'cell_volume_lambda': vol_data['growth_rate_hr'].values,
    'DNA_to_protein': dna_data['DNA_protein_ratio'].values,
    
    'lit_ms_lam': ms_data[ms_data['localization']=='envelope']['growth_rate_hr'].values,
    'lit_phi_rib': ms_data[ms_data['localization']=='ribosomal sector']['mass_frac'].values,
    'lit_phi_mem': ms_data[ms_data['localization']=='membrane']['mass_frac'].values,
    'lit_phi_cyt': ms_data[ms_data['localization']=='cytoplasm']['mass_frac'].values,
    'lit_phi_mem': ms_data[ms_data['localization']=='membrane']['mass_frac'].values,

    'ms_lam': data['growth_rate_hr'].values,
    'phi_rib': data['phi_rib'].values,
    'phi_cyt': data['phi_cyto'].values,
    'phi_mem': data['phi_mem'].values,
    'phi_peri': data['phi_peri'].values,
    'vol': data['volume'].values,
    'sa': data['surface_area'].values
}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

#%%
# Unpack the various samples, in an unfortunately complex manner. 
pars = ['rho_cyt', 'sigma_mem', 'density_ratio', 'predicted_SAV']
keys = ['date', 'strain', 'carbon_source', 'inducer_conc', 'replicate', 
        'growth_rate_hr', 'surface_to_volume']
pars_df = pd.DataFrame()
for i, p in enumerate(pars):
    post = samples.posterior[p].to_dataframe().reset_index()
    for g, d in post.groupby(f'{p}_dim_0'):
        _keys = data.iloc[g][keys]
        mean_val = d[p].mean()
        median_val = d[p].median()
        percs = np.percentile(d[p].values, [2.5, 97.5])
        _df = pd.DataFrame({'parameter': p,
                            'mean_val': mean_val,
                            'median_val': median_val,
                            'lower_95': percs[0],
                            'upper_95': percs[1]},
                            index=[0])
        for k, v in zip(keys, _keys):
            _df[k] = v
        pars_df = pd.concat([pars_df, _df])

#%%
# Unpack the literature data
ms_idx = ms_data[ms_data['localization']=='envelope'].copy().reset_index()
pars = ['lit_rho_cyt', 'lit_sigma_mem', 'lit_density_ratio']
keys = ['dataset_name', 'condition', 'growth_rate_hr'] 
lit_pars_df = pd.DataFrame()
for i, p in enumerate(pars):
    post = samples.posterior[p].to_dataframe().reset_index()
    for g, d in post.groupby(f'{p}_dim_0'):
        _keys = ms_idx.iloc[g][keys]
        mean_val = d[p].mean()
        median_val = d[p].median()
        percs = np.percentile(d[p].values, [2.5, 97.5])
        _df = pd.DataFrame({'parameter': p.split('lit_')[-1],
                            'mean_val': mean_val,
                            'median_val': median_val,
                            'lower_95': percs[0],
                            'upper_95': percs[1]},
                            index=[0])
        for k, v in zip(keys, _keys):
            _df[k] = v
        lit_pars_df = pd.concat([lit_pars_df, _df])

#%%
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
axes = {'ribosomal sector':ax[0,0], 
        'cytoplasm':ax[0, 1],
        'membrane': ax[0,2]}
for (g, d) in ms_data.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    for (par, _ax) in axes.items():
        _d = d[d['localization']==par]
        fmt = size.viz.style_point(g[0])
        fmt['ms'] = 4
        _ax.plot(_d['growth_rate_hr'], _d['mass_frac'], **fmt)

axes = {'rho_cyt': ax[1, 0],
        'sigma_mem': ax[1, 1],
        'density_ratio': ax[1, 2]}

for (g, d) in lit_pars_df.groupby(['dataset_name', 'condition', 'growth_rate_hr']):
    for (par, _ax) in axes.items():
        _d = d[d['parameter']==par]
        fmt = size.viz.style_point(g[0])
        fmt['ms'] = 4
        _ax.vlines(_d['growth_rate_hr'], _d['lower_95'], _d['upper_95'],
                   lw=0.5, color=cor['light_black'])
        _ax.plot(_d['growth_rate_hr'], _d['median_val'], **fmt)


wt_data = data[data['strain']=='wildtype']
axes = {'phi_rib': ax[0, 0],
        'phi_cyto': ax[0, 1],
        'phi_mem': ax[0, 2]}
fmt = {'marker':'o', 'markeredgewidth':1, 'markerfacecolor':'white',
       'markeredgecolor':cor['primary_black'], 'ms':4,
       'linestyle':'none'}
for (par, _ax) in axes.items():
    _ax.plot(wt_data['growth_rate_hr'], wt_data[par], **fmt)


concs = {'relA': {0:cor['pale_red'], 2:cor['primary_red'], 4:cor['red']},
         'meshI': {0:cor['pale_green'], 100:cor['primary_green']},
         'lacZ': {0:cor['pale_blue'], 1:cor['primary_blue'], 3:cor['blue'], 5:cor['dark_blue']}}
for g, d in data[data['strain']!='wildtype'].groupby(['strain', 'inducer_conc']):
    for (par, _ax) in axes.items():
        fmt['markerfacecolor'] = concs[g[0]][g[1]]
        _ax.vlines(d['growth_rate_hr'], d[par], d[par], lw=1, color=concs[g[0]][g[1]])
        _ax.plot(d['growth_rate_hr'], d[par], **fmt)

axes = {'rho_cyt': ax[1, 0],
        'sigma_mem': ax[1, 1],
        'density_ratio': ax[1, 2]}
fmt['markerfacecolor'] = 'white'
for (g, d) in pars_df[pars_df['strain']=='wildtype'].groupby('parameter'):
    if g == 'predicted_SAV':
        continue
    fmt['ms'] = 4
    axes[g].vlines(d['growth_rate_hr'], d['lower_95'],
                    d['upper_95'], lw=0.5, color=cor['primary_black'])
    axes[g].plot(d['growth_rate_hr'], d['mean_val'], **fmt)

for (g, d) in pars_df[(pars_df['strain'] != 'wildtype') & 
                      (pars_df['parameter']!='predicted_SAV')].groupby(['strain', 'inducer_conc', 'parameter']):
    fmt['markerfacecolor'] = concs[g[0]][g[1]]
    axes[g[-1]].vlines(d['growth_rate_hr'], d['lower_95'], d['upper_95'],lw=1, color=concs[g[0]][g[1]])
    axes[g[-1]].plot(d['growth_rate_hr'], d['median_val'], **fmt)

# Set the axes and labels
for a in ax.ravel():
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0,0].set_ylabel('ribosomal mass fraction', fontsize=6)
ax[0,1].set_ylabel('cytoplasmic mass fraction', fontsize=6)
ax[0,2].set_ylabel('membrane mass fraction', fontsize=6)

ax[1,0].set_ylabel('cytoplasmic density [fg/fL]', fontsize=6)
ax[1,1].set_ylabel('membrane density [fg/µm$^2$]', fontsize=6)
ax[1,2].set_ylabel('density ratio [µm$^{-1}$]', fontsize=6)

# Set the limits
ax[0,0].set_ylim([0, 0.35])
ax[0,1].set_ylim([0.5, 1])
ax[0,2].set_ylim([0, 0.25])
ax[1,0].set_ylim([0, 700])
ax[1,2].set_ylim([0, 250])
plt.savefig('./inferred_density_ratio.pdf', bbox_inches='tight')
#%%
# Plot predicted versus observed SAV
fig, ax = plt.subplots(1,1, figsize=(3,3))

_pred = pars_df[(pars_df['strain']=='wildtype') &
                (pars_df['parameter']=='predicted_SAV')]

fmt = {'marker':'o', 'markeredgewidth':1, 'markerfacecolor':'white',
       'markeredgecolor':cor['primary_black'], 'ms':4,
       'linestyle':'none'}
ax.vlines(_pred['surface_to_volume'], _pred['lower_95'], _pred['upper_95'],
            color=cor['primary_black'], lw=1)
ax.plot(_pred['surface_to_volume'], _pred['median_val'], **fmt)

# Plot the rels
concs = {'relA': {0:cor['pale_red'], 2:cor['primary_red'], 4:cor['red']},
         'meshI': {0:cor['pale_green'], 100:cor['primary_green']},
         'lacZ': {0:cor['pale_blue'], 1:cor['primary_blue'], 3:cor['blue'], 5:cor['dark_blue']}}
for (g, d) in pars_df[(pars_df['strain'] != 'wildtype') & 
                      (pars_df['parameter']=='predicted_SAV')].groupby(['strain', 'inducer_conc']):
    fmt['markerfacecolor'] = concs[g[0]][g[1]]
    ax.vlines(d['surface_to_volume'], d['lower_95'], d['upper_95'],lw=1, color=concs[g[0]][g[1]])
    ax.plot(d['surface_to_volume'], d['median_val'], **fmt)

ax.plot([4, 10], [4, 10], 'k--')
ax.set_xlabel('observed SAV [µm$^{-1}$]', fontsize=6)
ax.set_ylabel('predicted SAV [µm$^{-1}$]', fontsize=6)
plt.savefig('./predicted_vs_observed_SAV.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1,1, figsize=(3,3))
kappa = 150
data['predicted_SAV'] = kappa * data['phi_mem'] / (2 * (1 + 2.19 * data['phi_rib'] - data['phi_mem'] - data['phi_peri']))


fmt = {'marker':'o', 'markeredgewidth':1, 'markerfacecolor':'white',
       'markeredgecolor':cor['primary_black'], 'ms':4,
       'linestyle':'none'}
wt = data[data['strain']=='wildtype']
ax.plot(wt['surface_to_volume'], wt['predicted_SAV'], **fmt)

# Plot the rels
concs = {'relA': {0:cor['pale_red'], 2:cor['primary_red'], 4:cor['red']},
         'meshI': {0:cor['pale_green'], 100:cor['primary_green']},
         'lacZ': {0:cor['pale_blue'], 1:cor['primary_blue'], 3:cor['blue'], 5:cor['dark_blue']}}
for (g, d) in data[data['strain'] != 'wildtype'].groupby(['strain', 'inducer_conc']):
    fmt['markerfacecolor'] = concs[g[0]][g[1]]
    ax.plot(d['surface_to_volume'], d['predicted_SAV'], **fmt)

ax.plot([4, 10], [4, 10], 'k--')
ax.set_xlabel('observed SAV [µm$^{-1}$]', fontsize=6)
ax.set_ylabel('predicted SAV [µm$^{-1}$]', fontsize=6)
plt.savefig('./predicted_vs_observed_SAV_constant_kapp.pdf', bbox_inches='tight')



#%%
# Plot the empirics
fig, ax = plt.subplots(2, 2, figsize=(6, 3))
for g, d in prot_per_cell.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['ms'] = 4
    ax[0,0].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['ms'] = 4
    ax[0,1].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt)

for g, d in vol_data.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['ms'] = 4
    ax[1,0].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)

for g, d in dna_data.groupby('source'):
    fmt = size.viz.style_point(g)
    fmt['ms'] = 4
    ax[1,1].plot(d['growth_rate_hr'], d['DNA_protein_ratio'], **fmt)

# compute the summaries and store as a df
fit_df = pd.DataFrame()
pars = ['scan_vol', 'scan_sa', 'scan_prot', 'scan_theta']
for (i, p) in enumerate(pars):
    post = samples.posterior[p].to_dataframe().reset_index()
    for g, d in post.groupby(f'{p}_dim_0'):
        mean_val = d[p].mean()
        median_val = d[p].median()
        percs = np.percentile(d[p].values, [2.5, 97.5])
        _df = pd.DataFrame({'parameter': p,
                            'mean_val': mean_val,
                            'median_val': median_val,
                            'lower_95': percs[0],
                            'upper_95': percs[1],
                            'growth_rate_hr': lam_range[g]},
                            index=[0])
        fit_df = pd.concat([fit_df, _df])

axes = {'scan_prot': ax[0, 0],
        'scan_sa': ax[0, 1],
        'scan_vol': ax[1, 0],
        'scan_theta': ax[1, 1]}
for g, d in fit_df.groupby('parameter'):
    axes[g].fill_between(d['growth_rate_hr'], d['lower_95'], d['upper_95'],
                         color=cor['light_red'], alpha=0.5)
    axes[g].plot(d['growth_rate_hr'], d['median_val'], color=cor['red'])

for a in ax.ravel():
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

ax[0,0].set_ylabel('protein per cell [fg]', fontsize=6)
ax[0,0].set_yscale('log')
ax[0, 1].set_ylabel('surface area [µm$^2$]', fontsize=6)
ax[1, 0].set_ylabel('volume [µm$^3$]', fontsize=6)
ax[1, 0].set_yscale('log')
ax[1, 1].set_ylabel('DNA/protein ratio', fontsize=6)
plt.savefig('./parameter_fits.pdf', bbox_inches='tight')