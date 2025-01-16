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
model = cmdstanpy.CmdStanModel(stan_file='./fig2_inference.stan')

#%%
# Define the data dictionary for sampling
n_fit = 200
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
        'phi_cyto': ms_data['phi_cyto'].values,
        'phi_mem': ms_data['phi_mem'].values,
        'phi_peri': ms_data['phi_peri'].values,
        
        # Constants and fit variables
        'fit_lam': fit_lam,
        'W_PERI': 0.025
}

# Sample the model and extract to xarray. 
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

#%%
# Set a mapper for dimension to wildtype sample
mapper = {i:(ms_data.iloc[i]['carbon_source'], 
             ms_data.iloc[i]['replicate']) for i in range(len(ms_data))}

#%%
# Summarize the various quantities.
percs = [0.025, 0.975, 0.16, 0.84]  # Corresponding to 2 and 1 sigma (approx)


quants = ['tot_prot_per_cell', 'cyt_prot_per_cell',
          'cyt_rna_per_cell', 'peri_prot_per_cell',
          'mem_prot_per_cell', 'rho_cyt_prot',
          'rho_cyt_rna', 'rho_cyt_tot', 
          'rho_peri', 'sigma_mem']

# Group by 'carbon_source' and 'replicate' and compute the percentiles and mean for a specific column
def compute_group_stats(group, quantiles):
    quantiles_dict = {f'quantile_{q}': group.quantile(q) for q in quantiles}
    mean_dict = {'mean': group.mean()}
    return pd.Series({**mean_dict, **quantiles_dict})

# Apply the function to each group and each quantity
quant_df = []
for i, q in enumerate(quants):
    # Extract the posterior distribution
    post = samples.posterior[q].to_dataframe().reset_index()

    # Map conditions    
    post[['carbon_source', 'replicate']] = [mapper[v] for v in post[f'{q}_dim_0']]

    # Group and compute the percentiles and mean
    grouped = post.groupby(['carbon_source', 'replicate'])[q].apply(lambda x: compute_group_stats(x, percs)).unstack()

    # Flatten the hierarchical index and rename the columns
    grouped.columns = ['sig2_lower', 'sig2_upper', 'sig1_lower', 'sig1_upper', 'mean']
    grouped = grouped.reset_index()

    grouped['quantity'] = q
    grouped['strain'] = 'wildtype'
    grouped['replicate'] = grouped['replicate'].astype(int)
    quant_df.append(grouped)
quant_df = pd.concat(quant_df, sort=False)
quant_df.to_csv('../../data/mcmc/empirical_densities_summary.csv', index=False)

#%% 
# Summarize the posterior predictive checks
quants = ['prot_per_cell_ppc', 'rna_per_cell_ppc']
ppc_df = []
for i, q in enumerate(quants):
    post = samples.posterior[q].to_dataframe().reset_index()
    grouped = post.groupby(f'{q}_dim_0')[q].apply(
        lambda x: compute_group_stats(x, percs)).unstack().reset_index()
    grouped.drop(columns=[f'{q}_dim_0'], inplace=True)
    grouped.columns = ['sig2_lower', 'sig2_upper', 'sig1_lower', 'sig1_upper', 'mean']
    grouped['growth_rate_hr'] = fit_lam
    grouped['quantity'] = q
    ppc_df.append(grouped)
ppc_df = pd.concat(ppc_df, sort=False)
ppc_df.to_csv('../../data/mcmc/rna_protein_per_cell_ppc_summary.csv', index=False)
