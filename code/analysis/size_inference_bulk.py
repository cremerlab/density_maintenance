# %%
import arviz
import size.viz
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
data = data[(data['date'] != '2022-03-11_') &
            (data['date'] != '2022-09-29_')]


# Restrict to physical bounds
# data['width_median'] -= 0.3
# data = data[(data['width_median'] >= 0.1) & (data['width_median'] <= 2)]

hierarchical_model = cmdstanpy.CmdStanModel(
    stan_file='stan_models/hierarchical_size_inference_uncentered.stan')
simple_model = cmdstanpy.CmdStanModel(
    stan_file='stan_models/size_inference.stan')

# %%
# Define aspects of the hyper parameters to save
hyper_dfs = []
summary_dfs = []
hyper_vars = ['width_mu', 'length_mu', 'vol_mu', 'SAV_mu']
_hyper_vars = ['width_mu', 'length_mu', 'vol_mu', 'SAV_mu']
rename_cols = {'width_mu': 'width_um',
               'length_mu': 'length_um',
               'SAV_mu': 'SAV_inv_um',
               'vol_mu': 'volume_fL'}
inits = {'width_mu': 1,
         'length_mu': 1,
         'vol_mu': 1,
         'SAV_mu': 1}
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
    # , adapt_delta=0.9, iter_sampling=3500)
    samples = model.sample(data=data_dict)
    #    iter_warmup=1000,
    #    iter_sampling=5000)
    samples = arviz.from_cmdstanpy(samples)
    samples.to_netcdf(
        f'../../data/mcmc/wildtype_{g}_size_inference_object.netcdf')

    # Unpack hyperparameters
    hyper_df = samples.posterior[hyper_vars].to_dataframe().reset_index()
    lp = samples.sample_stats.lp
    mode_ind = np.where(lp == np.max(lp))

    # Compute various summary statistics of the hyper parameters
    for var in _hyper_vars:
        _percs = np.percentile(hyper_df[var], percs)
        _perc_df = pd.DataFrame([_percs], columns=perc_cols)
        _perc_df['mean'] = hyper_df[var].mean()
        _perc_df['median'] = hyper_df[var].median()
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
