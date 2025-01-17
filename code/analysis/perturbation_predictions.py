#%%
import numpy as np 
import pandas as pd 

# Load the experimental data and constrain to non-wildtype samples
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data = data[data['strain']!='wildtype']

# Load the parameter summaries for the density ratio from the mcmc. 
kappa = pd.read_csv('../../data/mcmc/theory_inference_parameter_summaries.csv')

# Define a function to compute the theory
def prediction(phi_mem: np.ndarray,
               phi_peri: np.ndarray, 
               phi_rib: np.ndarray, 
               kappa : float,
               beta: float = 1/0.4558) -> np.ndarray:
    numer = kappa * phi_mem
    denom = 2 * (1 + beta * phi_rib - phi_peri - phi_mem)
    return numer/denom

