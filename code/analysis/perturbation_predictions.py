#%%
import numpy as np 
import pandas as pd 

# Load the experimental data and constrain to non-wildtype samples
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')

# Load the parameter summaries for the density ratio from the mcmc. 
kappa = pd.read_csv('../../data/mcmc/theory_inference_parameter_summaries.csv')
kappa = kappa[kappa['quantity']=='kappa']

data = data[['replicate', 'strain', 'carbon_source', 'inducer_conc', 'surface_to_volume_inv_um', 'phi_mem', 'phi_peri', 'phi_rib', 'growth_rate_hr']]
data.rename(columns={'surface_to_volume_inv_um':'measured_SAV'}, inplace=True)

# Define a function to compute the theory
def prediction(kappa: float,
               beta: float = 2.19) -> np.ndarray:
    numer = kappa * data.phi_mem.values
    denom = 2 * (1 + beta * data.phi_rib.values - data.phi_peri.values - data.phi_mem.values)
    sav = numer / denom
    return sav

data['predicted_SAV_median'] = prediction(kappa['median'].values[0])
data['predicted_SAV_mean'] = prediction(kappa['mean'].values[0])
data['predicted_SAV_sig2_lower'] = prediction(kappa['sig2_lower'].values[0])
data['predicted_SAV_sig2_upper'] = prediction(kappa['sig2_upper'].values[0])
data['predicted_SAV_sig1_lower'] = prediction(kappa['sig1_lower'].values[0])
data['predicted_SAV_sig1_upper'] = prediction(kappa['sig1_upper'].values[0])
data.to_csv('../../data/mcmc/predicted_SAV_summary.csv', index=False)


