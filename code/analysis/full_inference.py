#%%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az
import size.viz
import matplotlib.pyplot as plt

# Compile the model
model = cmdstanpy.CmdStanModel(stan_file='./full_inference.stan')

#%%
# Load the literature size data
lit_size_data = pd.read_csv('../../data/literature/full_literature_size_data.csv')
lit_prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
lit_ms_data = pd.read_csv('../../data/collated/literature_mass_spec_aggregated.csv')
lit_ms_data = lit_ms_data[lit_ms_data['source'] != 'This Study']

# Merge the literature protein data and ours
prot_data = pd.read_csv('../../data/bulk_protein_per_cell.csv')
prot_data.rename(columns={'fg_prot_per_cell': 'fg_protein_per_cell',
                          'mean_growth_rate_hr': 'growth_rate_hr'}, inplace=True)
prot_data = prot_data[['growth_rate_hr', 'fg_protein_per_cell']]
prot_data['source'] = 'This Study'
prot_data = pd.concat([lit_prot_data, prot_data], sort=False)

# Load our coalesced data
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']
 
# Merge our size measurements and the literature
lit_size_data = lit_size_data[['source', 'growth_rate_hr', 'volume_um3', 'surface_area_um2']]
lit_size_data['replicate'] = 0
size_data = data[['growth_rate_hr', 'volume', 'surface_area', 'replicate']]
size_data['source'] = 'This Study'
size_data.rename(columns={'volume': 'volume_um3', 'surface_area': 'surface_area_um2'}, inplace=True)
size_data = pd.concat([lit_size_data, size_data], sort=False)


#%%
# Set the prediction range of growth rates and ribosome content
N_pred = 75
lam_range = np.linspace(0.005, 3, N_pred)
phi_rib_range = np.linspace(0, 0.45, N_pred)

# Define the data dictionary
data_dict = {
    'N_prot': len(prot_data),
    'N_size': len(size_data),
    'N_obs': len(data),
    'N_obs_lit': len(lit_ms_data),
    'N_pred': N_pred,

    'prot_growth_rate_hr': prot_data['growth_rate_hr'].values.astype(float),
    'prot_per_cell': prot_data['fg_protein_per_cell'].values.astype(float),
    'vol_growth_rate_hr': size_data['growth_rate_hr'].values.astype(float),
    'volume': size_data['volume_um3'].values.astype(float),
    'sa_growth_rate_hr': size_data['growth_rate_hr'].values.astype(float),
    'sa': size_data['surface_area_um2'].values.astype(float),

    'lit_growth_rate_hr': lit_ms_data['growth_rate_hr'].values.astype(float),
    'lit_phi_peri': lit_ms_data['phi_peri'].values.astype(float),
    'lit_phi_cyto': lit_ms_data['phi_cyt'].values.astype(float),
    'lit_phi_mem': lit_ms_data['phi_mem'].values.astype(float),
    'lit_phi_rib': lit_ms_data['phi_rib'].values.astype(float),

    'obs_growth_rate_hr': data['growth_rate_hr'].values.astype(float),
    'obs_sa': data['surface_area'].values.astype(float),
    'obs_sav': data['surface_to_volume'].values.astype(float),
    'obs_volume': data['volume'].values.astype(float),
    'obs_phi_peri': data['phi_peri'].values.astype(float),
    'obs_phi_cyto': data['phi_cyto'].values.astype(float),
    'obs_phi_mem': data['phi_mem'].values.astype(float),
    'obs_phi_rib': data['phi_rib'].values.astype(float),

    'DELTA_PERI': 0.025,
    'BETA_RIB': 1/0.4558,
    'pred_growth_rate_hr': lam_range,
    'pred_phi_rib': phi_rib_range
}

_samples = model.sample(data=data_dict, iter_sampling=4000)
samples = az.from_cmdstanpy(_samples)

#%%
sources = lit_ms_data['source'].values
growth_rates = lit_ms_data['growth_rate_hr'].values
condition = lit_ms_data['condition'].values

# Generate a summary dataframe from literature mass spec data
pars = ['m_peri_ppc', 'm_peri',
        'rho_cyto_ppc', 'rho_cyto',
        'rho_peri_ppc', 'rho_peri',
        'sigma_mem_ppc', 'sigma_mem'] 
quantity = ['m_peri_ppc', 'm_peri', 'rho_cyt_ppc', 'rho_cyt', 'rho_peri_ppc', 
        'rho_peri', 'sigma_mem_ppc', 'sigma_mem']
ms_densities = pd.DataFrame([])
for (p, q) in zip(pars, quantity): 
    _d = samples.posterior[f'lit_{p}'].to_dataframe().reset_index()
    for i in range(len(sources)):
        _d.loc[_d[f'lit_{p}_dim_0']==i, 'source'] = sources[i]
        _d.loc[_d[f'lit_{p}_dim_0']==i, 'growth_rate_hr'] = growth_rates[i]
        _d.loc[_d[f'lit_{p}_dim_0']==i, 'condition'] = condition[i]
    _d.rename(columns={f'lit_{p}': q}, inplace=True)
    percs = size.viz.compute_percentiles(_d, q, ['source', 'growth_rate_hr', 'condition'])
    percs['replicate'] = 0
    ms_densities = pd.concat([ms_densities, percs])

# %%
# Generate a summary dataframe from our mass spec data
growth_rates = data['growth_rate_hr'].values
ondition = data['carbon_source'].values
replicate = data['replicate'].values
for (p, q) in zip(pars, quantity): 
    _d = samples.posterior[p].to_dataframe().reset_index()
    for i in range(len(growth_rates)):
        _d.loc[_d[f'{p}_dim_0']==i, 'growth_rate_hr'] = growth_rates[i]
        _d.loc[_d[f'{p}_dim_0']==i, 'condition'] = condition[i]
        _d.loc[_d[f'{p}_dim_0']==i, 'replicate'] = replicate[i]
    _d.rename(columns={p: q}, inplace=True)
    percs = size.viz.compute_percentiles(_d, q, ['growth_rate_hr', 'condition', 'replicate'])
    percs['source'] = 'This Study'
    ms_densities = pd.concat([ms_densities, percs])

ms_densities.to_csv('./output/mass_spec_densities_summary.csv', index=False)

#%%
# Generate a wide dataframe with measured/predicted values for each quantity
pars = ['pred_sav_ppc', 'pred_sav']
parnames = ['sav_ppc', 'sav']
carbon_sources = data['carbon_source'].values
growth_rates = data['growth_rate_hr'].values
replicates = data['replicate'].values
pred_sav = pd.DataFrame([])
for (p, n) in zip(pars, parnames):
    _d = samples.posterior[p].to_dataframe().reset_index()
    for i in range(len(carbon_sources)):
        _d.loc[_d[f'{p}_dim_0']==i, 'condition'] = carbon_sources[i]
        _d.loc[_d[f'{p}_dim_0']==i, 'growth_rate_hr'] = growth_rates[i]
        _d.loc[_d[f'{p}_dim_0']==i, 'replicate'] = replicates[i]
    _d.rename(columns={p: n}, inplace=True)
    percs = size.viz.compute_percentiles(_d, n, ['growth_rate_hr', 'condition', 'replicate'])
    percs['source'] = 'This Study'
    pred_sav = pd.concat([pred_sav, percs])
pred_sav.to_csv('./output/predicted_sav_summary.csv', index=False)

#%%
# Generate a summary dataframe of the fits
pars = ['fit_prot', 'fit_volume', 'fit_sa']
quantity = ['protein_per_cell', 'volume_um3', 'surface_area_um2']
fits = pd.DataFrame([])
for (p, q) in zip(pars, quantity):
    _d = samples.posterior[p].to_dataframe().reset_index()
    for i, ell in enumerate(lam_range):
        _d.loc[_d[f'{p}_dim_0']==i, 'growth_rate_hr'] = ell
    _d.rename(columns={p: q}, inplace=True)
    percs = size.viz.compute_percentiles(_d, q, ['growth_rate_hr'])
    fits = pd.concat([fits, percs])
fits.to_csv('./output/size_relation_fits_summary.csv', index=False)

# Generate a dataframe of the phi_rib scaling fits and SAV predictions
pars = ['phi_mem_ppc', 'phi_mem', 'phi_peri_ppc', 'phi_peri', 'sav_ppc', 'sav']
preds = pd.DataFrame([])
for p in pars:
    if 'sav' in p:
        prefix = 'theory_'
    else:
        prefix='fit_'
    _d = samples.posterior[f'{prefix}{p}'].to_dataframe().reset_index()
    for i, phi_rib in enumerate(phi_rib_range):
        _d.loc[_d[f'{prefix}{p}_dim_0']==i, 'phi_rib'] = phi_rib
    _d.rename(columns={f'{prefix}{p}': p}, inplace=True)
    percs = size.viz.compute_percentiles(_d, p, ['phi_rib'])
    preds = pd.concat([preds, percs])
preds.to_csv('./output/phi_rib_scaling_fits_summary.csv', index=False)

#%%
# Scaling 
pars = ['beta_0_phi_mem', 'beta_1_phi_mem', 'phi_mem_sigma',
        'beta_0_phi_peri', 'beta_1_phi_peri', 'phi_peri_sigma']
dist = samples.posterior[pars].to_dataframe().reset_index()
dist.to_csv('./output/phi_rib_scaling_fits_samples.csv', index=False)

#%%
# Save full distribution of mean densites
pars = ['kappa', 'sav_sigma']
dist = samples.posterior[pars].to_dataframe().reset_index()
dist = dist[pars]
dist.rename(columns=mapper, inplace=True)
dist.to_csv('./output/kappa_density_samples.csv', index=False)
#%%
cor, pal = size.viz.matplotlib_style()
# Generate a diagnostic plot of the fits
fig, ax = plt.subplots(1, 3, figsize=(6, 2))

for g, d in lit_size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)
    ax[1].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt)

for g, d in lit_prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

# Plot the medians
axes = {'volume_um3':ax[0], 'surface_area_um2':ax[1], 'protein_per_cell':ax[2]}
for (q, a) in axes.items():
    _fits = fits[fits['quantity']==q]
    meds = _fits[_fits['interval']=='median']
    perc = _fits[_fits['interval']=='95%']
    a.fill_between(perc['growth_rate_hr'], perc['lower'], perc['upper'], color=cor['primary_black'], alpha=0.2)
    a.plot(meds['growth_rate_hr'], meds['lower'], color=cor['primary_black'], alpha=1)

