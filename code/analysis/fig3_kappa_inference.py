#%%
import numpy as np 
import pandas as pd
import cmdstanpy 
import arviz as az

# Load the dataset
data = pd.read_csv('../../data/compiled_measurements.csv')
data = data[data['strain']=='wildtype']

# Load and compile the inferential model
model = cmdstanpy.CmdStanModel(stan_file='fig3_kappa_inference.stan')

#%%
# Define the range of ribosomal allocation over which to draw the ppcs
phi_rib_range = np.linspace(0, 0.45, 100)
# Define the data dictionary
data_dict = {

    # Define the dimensional information
    'N_obs': len(data),
    'N_fit': len(phi_rib_range),

    # Define the observed data
    'obs_phi_rib': data['phi_rib'].values,
    'obs_phi_mem': data['phi_mem'].values,
    'obs_phi_peri': data['phi_peri'].values,
    'obs_sav': data['sav_inv_um'].values,

    # Define the ribosomal allocation range
    'phi_rib_range': phi_rib_range,

    # Define the conversion factor for ribosomal allocation to RNA-to-protein
    'BETA': 1/0.4558,
}

# Sample the model and extract
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

#%%
# Extract the estimates for the parameters. 
pars = ['kappa', 'sigma', 'beta_0_phi_mem', 'beta_1_phi_mem', 'sigma_phi_mem',
        'beta_0_phi_peri', 'beta_1_phi_peri', 'sigma_phi_peri']
post = samples.posterior[pars].to_dataframe().reset_index()
melted = post.melt(value_vars=pars)

# Generate a tidy-longform with mean and percentiles
percs = [(2.5, 97.5), (18, 62)]
labels = ['2sig_', '1sig_']
dfs = []
for g, d in melted.groupby('variable'):
    _df = {'quantity': g,
           'median_value': d['value'].median()}
    for p, ell in zip(percs, labels):
        lower, upper = np.percentile(d['value'].values, p) 
        _df[f'{ell}lower'] = lower
        _df[f'{ell}upper'] = upper
    dfs.append(pd.DataFrame(_df, index=[0]))
param_df = pd.concat(dfs, sort=False)
param_df.to_csv('../../data/mcmc/fig3_parameter_summaries.csv', index=False)

#%%
# Extract the ppcs for the theorya nd the trend fits
quants = ['phi_mem_mu', 'phi_mem_ppc', 'phi_peri_mu', 'phi_peri_ppc', 
        'theory_mu', 'theory_ppc']

# Compute the percentiles for each parameter
dfs = []
for q in quants:
    post = samples.posterior[[q, f'{q}_dim_0']].to_dataframe().reset_index()
    for g, d in post.groupby(f'{q}_dim_0'):
        _df = {'quantity': q,
               'phi_rib': phi_rib_range[g],
               'median_val': np.percentile(d[q].values, [50, 50])[0]}
        for p, ell in zip(percs, labels):
            lower, upper = np.percentile(d[q].values, p)
            _df[f'{ell}lower'] = lower
            _df[f'{ell}upper'] = upper
        dfs.append(pd.DataFrame(_df, index=[0]))
fit_df = pd.concat(dfs, sort=False)
fit_df.to_csv('../../data/mcmc/fig3_fits.csv', index=False)