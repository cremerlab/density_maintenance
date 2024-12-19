#%%
import numpy as np 
import pandas as pd 
import cmdstanpy
import arviz as az

# Load the model
model = cmdstanpy.CmdStanModel(stan_file='fig2_mass_density_inference.stan')

#%%
# Load our measurements
data = pd.read_csv('../../data/compiled_measurements.csv')
data = data[data['strain']=='wildtype']
data['source'] = 'This Study'

# Load literature mass spec measurements
lit_ms_data = pd.read_csv('../../data/collated/compiled_literature_allocation_assignments_wide.csv')

# Load the literature size data
lit_size_data = pd.read_csv('../../data/literature/full_literature_size_data.csv')

# Load the protein per cell data and combine for inference
lit_prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
prot_data = pd.read_csv('../../data/bulk_protein_per_cell.csv')

# Concatenate protein per cell measurements for inference
prot_per_cell = np.concatenate([lit_prot_data['fg_protein_per_cell'].values,
                                   prot_data['fg_prot_per_cell'].values])
prot_per_cell_lam = np.concatenate([lit_prot_data['growth_rate_hr'].values,
                                    prot_data['mean_growth_rate_hr'].values])

# Concatenate size measurements for inference
volume = np.concatenate([lit_size_data['volume_um3'].values,
                         data['volume_fL'].values])
surface_area = np.concatenate([lit_size_data['surface_area_um2'].values,
                               data['surface_area_um2'].values])
size_lam = np.concatenate([lit_size_data['growth_rate_hr'].values,
                           data['growth_rate_hr'].values])

# Concatenate ms allocation measurements for empirical calculation
ms_lam = np.concatenate([lit_ms_data['growth_rate_hr'].values,
                         data['growth_rate_hr'].values])
ms_phi_cyto = np.concatenate([lit_ms_data['phi_cyto'].values,
                              data['phi_cyto'].values])
ms_phi_peri = np.concatenate([lit_ms_data['phi_peri'].values,
                              data['phi_peri'].values])
ms_phi_mem = np.concatenate([lit_ms_data['phi_mem'].values,
                             data['phi_mem'].values])
ms_phi_rib = np.concatenate([lit_ms_data['phi_rib'].values, data['phi_rib'].values])
ms_source = np.concatenate([lit_ms_data['source'].values,
                            data['source'].values])
# Set the details for the fit
fit_lam = np.linspace(0, 2.5, 100)
#%%
# Set the data dictionary
data_dict = {
    # Dimensionality
   'N_prot': len(prot_per_cell),
   'N_size': len(volume),
   'N_ms': len(ms_lam),
   'N_fit': len(fit_lam),
   'N_obs': len(data),

   # Independent data
   'prot_per_cell': prot_per_cell,
   'surface_area': surface_area,
   'volume': volume,

   # Dependent data
   'prot_per_cell_lam': prot_per_cell_lam,
   'size_lam': size_lam,

   # Our experimental measurements for easy calculation of densities
   'obs_volume': data['volume_fL'].values,
   'obs_surface_area': data['surface_area_um2'].values,
   'obs_phi_cyto': data['phi_cyto'].values,
   'obs_phi_mem': data['phi_mem'].values,
   'obs_phi_peri': data['phi_peri'].values,
   'obs_phi_rib': data['phi_rib'].values,
   'obs_lam': data['growth_rate_hr'].values,

   # Generative modeling
   'phi_cyto': ms_phi_cyto,
   'phi_mem': ms_phi_mem,
   'phi_peri': ms_phi_peri,
   'phi_rib': ms_phi_rib, 
   'ms_lam': ms_lam,
   'fit_lam': fit_lam,
   
   # Constant for periplasm width
   'W_PERI': 0.025,
   'BETA': 1/0.4558,
   'rRNA_FRAC': 0.8
}

# Sample the model 
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

#%%
# Compute the mean and percentiles for each mass spec point. 
quantities = ['emp_M_cyto', 'emp_M_peri', 'emp_M_mem', 'emp_rho_cyto', 
              'emp_rho_peri', 'emp_sigma_mem', 'emp_rho_rib', 
              'emp_rho_rrna', 'emp_rho_rna', 'emp_rho_biomass']
dfs = [] 
for q in quantities:
    post = samples.posterior[[q, f'{q}_dim_0']].to_dataframe().reset_index()
    for g, d in post.groupby([f'{q}_dim_0']):
        mean_val = np.mean(d[q])
        percs = np.percentile(d[q].values, (2.5, 97.5))
        dfs.append(pd.DataFrame({'quantity':q[4:],
                           'source':ms_source[g],
                           'growth_rate_hr': ms_lam[g],
                           'mean_val':mean_val,
                           'lower':percs[0],
                           'upper':percs[1]},
                           index=[0])) 
mass_densities = pd.concat(dfs, sort=False)
mass_densities = mass_densities[mass_densities['source'] != 'This Study']
# Compute the mean and percentiles for our experimental measurements
quantities = ['obs_M_cyto', 'obs_M_peri', 'obs_M_mem', 'obs_rho_cyto', 
              'obs_rho_peri', 'obs_sigma_mem', 'obs_rho_rib', 'obs_rho_rrna',
              'obs_rho_rna', 'obs_rho_biomass']
dfs = []
for q in quantities:
    post = samples.posterior[[q, f'{q}_dim_0']].to_dataframe().reset_index()
    for g, d in post.groupby([f'{q}_dim_0']):
         mean_val = np.mean(d[q])
         percs = np.percentile(d[q].values, (2.5, 97.5))
         dfs.append(pd.DataFrame({'quantity':q[4:],
                           'source':'This Study',
                           'growth_rate_hr':data['growth_rate_hr'].values[g],
                           'mean_val':mean_val,
                           'lower':percs[0],
                           'upper':percs[1]},
                           index=[0])) 
_mass_densities = pd.concat(dfs, sort=False)
mass_densities = pd.concat([mass_densities, _mass_densities], sort=False)
mass_densities.to_csv('../../data/mcmc/fig2_empirical_quantities_summary.csv',
                      index=False)

#%%
# Compute the best fit over the growth relations with 95% percentiles
quantities = ['prot_per_cell_ppc', 'volume_ppc', 'surface_area_ppc']
dfs = []
for q in quantities:
    post = samples.posterior[[q, f'{q}_dim_0']].to_dataframe().reset_index()
    for g, d in post.groupby(f'{q}_dim_0'):
        mean_val = np.mean(d[q])
        percs = np.percentile(d[q], [2.5, 97.5])
        dfs.append(pd.DataFrame({'quantity':q,
                                 'growth_rate_hr': fit_lam[g],
                                 'mean_val':mean_val,
                                 'lower':percs[0],
                                 'upper':percs[1]},
                                 index=[0]))

empirical_fits = pd.concat(dfs, sort=False)
empirical_fits.to_csv('../../data/mcmc/fig2_ppc_summary.csv', index=False)


    
