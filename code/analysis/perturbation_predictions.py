#%%
import numpy as np 
import pandas as pd 

# Load the experimental data and constrain to non-wildtype samples
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data = data[data['strain']!='wildtype']

# Load the parameter summaries for the density ratio from the mcmc. 
kappa = pd.read_csv('../../data/mcmc/theory_inference_parameter_summaries.csv')
kappa = kappa[kappa['quantity']=='kappa']



data = data[['replicate', 'strain', 'carbon_source', 'inducer_conc', 'surface_to_volume_inv_um', 'phi_mem', 'phi_peri', 'phi_rib']]
data.rename(columns={'surface_to_volume_inv_um':'measured_SAV'}, inplace=True)

# Define a function to compute the theory
def prediction(kappa: float,
               beta: float = 1/0.4558) -> np.ndarray:
    numer = kappa * data.phi_mem.values
    denom = 2 * (1 + beta * data.phi_rib.values - data.phi_peri.values - data.phi_mem.values)
    sav = numer / denom
    return sav


data['predicted_SAV_mean'] = prediction(kappa['mean'].values[0])
data['predicted_SAV_2sig_lower'] = prediction(kappa['2sig_lower'].values[0])
data['predicted_SAV_2sig_upper'] = prediction(kappa['2sig_upper'].values[0])
data['predicted_SAV_1sig_lower'] = prediction(kappa['1sig_lower'].values[0])
data['predicted_SAV_1sig_upper'] = prediction(kappa['1sig_upper'].values[0])
data.to_csv('../../data/mcmc/perturbation_predicted_SAV_summary.csv', index=False)


