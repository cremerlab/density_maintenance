#%%
import pandas as pd
import numpy as np 
import cmdstanpy 
import arviz as az

# Load the data and restrict to wildtype as necessary
ms_data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
ms_data = ms_data[ms_data['strain']=='wildtype']
rp_data = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')

# Load and compile the inference model
model = cmdstanpy.CmdStanModel(stan_file='./empirical_density_inference.stan')

#%%
# Define the data dictionary for sampling
n_fit = 30
fit_lam = np.linspace(0, 3, n_fit)
data_dict = {
        # Dimensional information
        'N_rp': len(rp_data),
        'N_obs': len(ms_data),
        'N_fit': n_fit,

        # Independent variables
        'prot_per_cell': rp_data['fg_protein_per_cell'].values,
        'rna_per_cell': rp_data['fg_rna_per_cell'].values,
        'surface_area': ms_data['surface_area_um2'].values,
        'volume': ms_data['volume_fL'].values,
        
        # Dependent variables
        'rp_lam': rp_data['growth_rate_hr'].values,
        'lam': ms_data['growth_rate_hr'].values,

        # Allocation parameters for empricial quantity calculation 
        'psi_cyto': ms_data['psi_cyto'].values,
        'psi_mem': ms_data['psi_mem'].values,
        'psi_peri': ms_data['psi_peri'].values,
        
        # Constants and fit variables
        'fit_lam': fit_lam,
        'W_PERI': 0.025 # in um
}

# Sample the model and extract to xarray. 
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

#%%
# Set a mapper for dimension to wildtype sample
mapper = {i:(ms_data.iloc[i]['carbon_source'], 
             ms_data.iloc[i]['replicate'], 
             ms_data.iloc[i]['growth_rate_hr'],
             ms_data.iloc[i]['strain'],
             ms_data.iloc[i]['inducer_conc']) for i in range(len(ms_data))}

#%%
# Summarize the various quantities.
percs = [2.5, 97.5, 16, 84] # Corresponds to lower and upper bounds of approx 2 and 1 sigma
quants = ['tot_prot_per_cell', 'cyt_prot_per_cell',
          'cyt_rna_per_cell', 'cyt_tot_per_cell','peri_prot_per_cell',
          'mem_prot_per_cell', 'rho_cyt_prot',
          'rho_cyt_rna', 'rho_cyt_tot', 
          'rho_peri', 'sigma_mem', 'empirical_kappa']

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
quant_df.to_csv('../../data/mcmc/empirical_densities_summary.csv', index=False)

#%% 
# Summarize the posterior predictive checks
quants = ['prot_per_cell_ppc', 'rna_per_cell_ppc']
ppc_df = []
dfs = []
for i, q in enumerate(quants):
    post = samples.posterior[q].to_dataframe().reset_index()
    for g, d in post.groupby(f'{q}_dim_0'):
        quants = np.percentile(d[q], percs)
        _df = pd.DataFrame({'quantity': q,
                            'growth_rate_hr': fit_lam[g],
                            'mean': d[q].mean(),
                            'median': d[q].median(),
                            'sig2_lower': quants[0],
                            'sig2_upper': quants[1],
                            'sig1_lower': quants[2],
                            'sig1_upper': quants[3]},
                            index=[0])
        dfs.append(_df)
ppc_df = pd.concat(dfs, sort=False)
ppc_df.to_csv('../../data/mcmc/rna_protein_per_cell_ppc_summary.csv', index=False)
