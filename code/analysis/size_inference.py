# %%
import arviz
import size.viz
import matplotlib.pyplot as plt
import cmdstanpy
import pandas as pd
import numpy as np
import tqdm
# %%
cor, pal = size.viz.matplotlib_style()
PERIPLASMIC_DIAM = 0.025  # in microns

# Load data and compile the stan model
data = pd.read_csv(
    '../processing/microscopy/wildtype_size_measurement/output/wildtype_size_measurements.csv')
hierarchical_model = cmdstanpy.CmdStanModel(
    stan_file='stan_models/hierarchical_size_inference.stan')
simple_model = cmdstanpy.CmdStanModel(
    stan_file='stan_models/size_inference.stan')

# %%

# Define aspects of the hyper parameters to save
hyper_dfs = []
summary_dfs = []
hyper_vars = ['width_mu', 'length_mu', 'vol_mu', 'sa_mu', 'ar_mu']
rename_cols = {'width_mu': 'width_um',
               'length_mu': 'length_um',
               'vol_mu': 'volume_fL',
               'sa_mu': 'surface_area_um2',
               'ar_mu': 'aspect_ratio'}

# Define the percentiles to compute for each hyperparameter
percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']

# Iterate through each carbon source and perform the inference.
for g, d in tqdm.tqdm(data.groupby(['carbon_source'])):
    # Assign replicate identifiers
    d = d.copy()
    d['idx'] = d.groupby(['date']).ngroup() + 1

    # Assemble the data dictionary
    data_dict = {
        'J': d['idx'].max(),
        'N': len(d),
        'idx': d['idx'].values.astype(int),
        'widths': d['width_median'].values.astype(float),
        'lengths': d['length'].values.astype(float),
        'periplasmic_diam': PERIPLASMIC_DIAM}

    # Sample the model and save full output to disk
    if data_dict['J'] > 1:
        model = hierarchical_model
    else:
        model = simple_model
    samples = model.sample(data=data_dict, iter_warmup=5000, iter_sampling=2000,
                           adapt_delta=0.9)
    samples = arviz.from_cmdstanpy(samples)
    samples.to_netcdf(f'../../data/mcmc/wildtype_{g}_size_inference_object.netcdf')

    # Unpack hyperparameters
    hyper_df = samples.posterior[hyper_vars].to_dataframe().reset_index()
    lp = samples.sample_stats.lp
    mode_ind = np.where(lp == np.max(lp))

    # Compute various summary statistics of the hyper parameters
    for var in hyper_vars:
        _percs = np.percentile(hyper_df[var], percs)
        _perc_df = pd.DataFrame([_percs], columns=perc_cols)
        _perc_df['mean'] = hyper_df[var].mean()
        _perc_df['median'] = hyper_df[var].median()
        # _perc_df['mode'] = float(samples.posterior[var][mode_ind].values[0])
        _perc_df['parameter'] = rename_cols[var]
        _perc_df['carbon_source'] = g
        _perc_df['strain'] = 'wildtype'
        summary_dfs.append(_perc_df)

    # Format and save dataframes
    hyper_df['carbon_source'] = g
    hyper_df['strain'] = 'wildtype'
    hyper_df.rename(columns=rename_cols, inplace=True)
    hyper_df.drop(columns=['chain', 'draw'], inplace=True)
    hyper_dfs.append(hyper_df)

# Concatenate and save the hyperparameters dataframes
hyper_df = pd.concat(hyper_dfs, sort=False)
summary_df = pd.concat(summary_dfs, sort=False)

# %%
hyper_df.to_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_samples.csv', index=False)
summary_df.to_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv', index=False)

# %%
