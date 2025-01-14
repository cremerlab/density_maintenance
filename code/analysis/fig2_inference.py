#%%
import pandas as pd
import numpy as np 
import cmdstanpy 
import arviz as az

data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
model = cmdstanpy.CmdStanModel(stan_file='./fig2_inference.stan')