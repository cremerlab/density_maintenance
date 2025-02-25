#%%
import numpy as np 
import pandas as pd 
import cmdstanpy
import arviz as az
model = cmdstanpy.CmdStanModel(stan_file='./perturbation_density_inference.stan')

#%%
# Load the relevant datasets 
prot = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
wildtype = data[data['strain'] == 'wildtype']
perturb = data[data['strain']!= 'wildtype']

# Set the growth rate rate
lam_range = np.linspace(0, 2.5, 50)

# Define the data dictionary
data_dict = {
    'N_prot': len(prot),
    'N_size': len(wildtype),
    'N_ppc': len(lam_range),
    'N_obs': len(perturb),

    'protein_per_cell': prot['fg_protein_per_cell'].values,
    'protein_per_cell_lam': prot['growth_rate_hr'].values,
    'volume': wildtype['volume_fL'].values,
    'volume_lam': wildtype['growth_rate_hr'].values,

    'phi_rib': perturb['phi_rib'].values,
    'phi_cyto': perturb['phi_cyto'].values,
    'phi_peri': perturb['phi_peri'].values,
    'phi_mem': perturb['phi_mem'].values,
    'obs_surface_area': perturb['surface_area_um2'].values,
    'obs_volume': perturb['volume_fL'].values,

    'ppc_lam_range': lam_range,
    'BETA_RIB': 1/0.4558,
    'W_PERI': 0.025
}

_samples = model.sample(data_dict)
samples = az.from_cmdstanpy(_samples)

#%%
# Set a mapper for dimension to perturbation sample
mapper = {i:(perturb.iloc[i]['carbon_source'], 
             perturb.iloc[i]['replicate'], 
             perturb.iloc[i]['growth_rate_hr'],
             perturb.iloc[i]['strain'],
             perturb.iloc[i]['inducer_conc']) for i in range(len(perturb))}

# Summarize the various quantities.
percs = [2.5, 97.5, 16, 84] # Corresponds to lower and upper bounds of approx 2 and 1 sigma
quants = ['prot_per_cell', 'rna_per_cell', 'cyt_tot_per_cell',
          'peri_prot_per_cell', 'mem_prot_per_cell',
          'rho_cyt', 'rho_peri', 'sigma_mem', 'empirical_kappa']

# Apply the function to each group and each quantity
quant_df = []
for i, q in enumerate(quants):
    # Extract the posterior distribution
    post = samples.posterior[q].to_dataframe().reset_index()
    for g, d in post.groupby(f'{q}_dim_0'):
        quants = np.percentile(d[q], percs)
        _df = pd.DataFrame({'quantity': q,
                            'strain': mapper[g][3],
                            'inducer_conc': mapper[g][4],
                            'carbon_source': mapper[g][0],
                            'replicate': mapper[g][1],
                            'growth_rate_hr': mapper[g][2],
                            'mean': d[q].mean(),
                            'median': d[q].median(),
                            'sig2_lower': quants[0],
                            'sig2_upper': quants[1],
                            'sig1_lower': quants[2],
                            'sig1_upper': quants[3]},
                            index=[0])
        quant_df.append(_df)
quant_df = pd.concat(quant_df, sort=False)
quant_df.to_csv('../../data/mcmc/perturbation_empirical_densities_summary.csv', index=False)


#%%
# Summarize the posterior predictive checks
quants = ['volume_ppc', 'protein_per_cell_ppc']
ppc_df = []
dfs = []
for i, q in enumerate(quants):
    post = samples.posterior[q].to_dataframe().reset_index()
    for g, d in post.groupby(f'{q}_dim_0'):
        quants = np.percentile(d[q], percs)
        _df = pd.DataFrame({'quantity': q,
                            'growth_rate_hr': lam_range[g],
                            'mean': d[q].mean(),
                            'median': d[q].median(),
                            'sig2_lower': quants[0],
                            'sig2_upper': quants[1],
                            'sig1_lower': quants[2],
                            'sig1_upper': quants[3]},
                            index=[0])
        dfs.append(_df)
ppc_df = pd.concat(dfs, sort=False)
ppc_df.to_csv('../../data/mcmc/volume_protein_per_cell_ppc_summary.csv', index=False)

#%%
#Save the samples of the estimated constant protein per cell 
df = samples.posterior.rho_prot_mu.to_dataframe().reset_index()[['rho_prot_mu']]
df.to_csv('../../data/mcmc/estimated_constant_total_protein_density.csv', index=False)