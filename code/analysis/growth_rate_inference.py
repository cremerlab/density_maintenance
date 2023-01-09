# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import size.viz
import arviz
import tqdm
colors, palette = size.viz.matplotlib_style()

# Load the processed data
data = pd.read_csv(
    '../../data/growth_curves/wt_growth_measurements_processed.csv')

# Load and compile the inference model
hier_model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/hierarchical_growth_inference.stan')
# Load and compile the inference model
simple_model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/growth_inference.stan')

# %%
lam_df = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby(['carbon_source'])):
    # Assign grouping for hierarchical analysis
    d = d.copy()
    d['J_idx'] = (d.groupby(
        'date').ngroup() + 1).astype(int)

    # Define the data dictionary
    data_dict = {'J': d['J_idx'].max(),
                 'N': len(d),
                 'idx': d['J_idx'].values,
                 'elapsed_time': d['elapsed_time_hr'].values.astype(float),
                 'optical_density': d['od_600nm'].values.astype(float)}
    if data_dict['J'] > 1:
        model = hier_model
    else:
        model = simple_model
    samples = model.sample(data=data_dict, adapt_delta=0.95)

    samples = arviz.from_cmdstanpy(samples)
    samples_df = samples.posterior.to_dataframe().reset_index()
    print('Processing growth rate....')
    lam = [d['mu'].values[0] for _, d in samples_df.groupby(['chain', 'draw'])]
    print('done!')
    _df = pd.DataFrame([])
    _df['growth_rate_hr'] = lam
    _df['carbon_source'] = g
    lam_df = pd.concat([lam_df, _df], sort=False)

# %%

# Compute the percentiles of the samples.
5percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']
hyper_vars = ['mu', 'mu_1']
summary_df = pd.DataFrame([])
for g, d in lam_df.groupby(['carbon_source']):
    _percs = np.percentile(d['growth_rate_hr'], percs)
    _perc_df = pd.DataFrame([_percs], columns=perc_cols)
    _perc_df['mean'] = d['growth_rate_hr'].mean()
    _perc_df['median'] = d['growth_rate_hr'].median()
    _perc_df['carbon_source'] = g
    _perc_df['strain'] = 'wildtype'
    summary_df = pd.concat([summary_df, _perc_df])
summary_df.to_csv(
    '../../data/mcmc/wildtype_growth_rate_summary.csv', index=False)
