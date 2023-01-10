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
    '../processing/microscopy/size_measurement/output/compiled_size_measurements.csv')
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
hyper_vars = ['width_mu', 'length_alpha', 'length_beta', 'vol_mu', 'SAV_mu']
_hyper_vars = ['width_mu', 'length_mu', 'vol_mu', 'SAV_mu']
rename_cols = {'width_mu': 'width_um',
               'length_mu': 'length_um',
               'length_alpha': 'length_alpha',
               'length_beta': 'length_beta',
               'SAV_mu': 'SAV_inv_um',
               'vol_mu': 'volume_fL'}

# Define the percentiles to compute for each hyperparameter
percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']

# Iterate through each carbon source and perform the inference.
for g, d in tqdm.tqdm(data.groupby(['carbon_source', 'strain', 'inducer', 'inducer_conc', 'temperature_C'])):
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
    # if data_dict['J'] > 1:
    model = hierarchical_model
    # else:
    # model = simple_model
    samples = model.sample(data=data_dict, adapt_delta=0.95,
                           iter_warmup=1000, iter_sampling=2000)

    samples = arviz.from_cmdstanpy(samples)

    # Unpack hyperparameters
    hyper_df = samples.posterior[hyper_vars].to_dataframe().reset_index()
    hyper_df['length_mu'] = hyper_df['length_alpha'].values / \
        hyper_df['length_alpha'].values
    lp = samples.sample_stats.lp
    mode_ind = np.where(lp == np.max(lp))

    # Compute various summary statistics of the hyper parameters
    for var in _hyper_vars:
        _percs = np.percentile(hyper_df[var], percs)
        _perc_df = pd.DataFrame([_percs], columns=perc_cols)
        _perc_df['mean'] = hyper_df[var].mean()
        _perc_df['median'] = hyper_df[var].median()
        _perc_df['parameter'] = rename_cols[var]
        _perc_df['carbon_source'] = g[0]
        _perc_df['strain'] = g[1]
        _perc_df['inducer'] = g[2]
        _perc_df['inducer_conc'] = g[3]
        _perc_df['temperature_C'] = g[4]
        summary_dfs.append(_perc_df)

    # Format and save dataframes
    hyper_df['carbon_source'] = g[0]
    hyper_df['strain'] = g[1]
    hyper_df['inducer'] = g[2]
    hyper_df['inducer_conc'] = g[3]
    hyper_df['temperature_C'] = g[4]
    hyper_df.rename(columns=rename_cols, inplace=True)
    hyper_df.drop(columns=['chain', 'draw'], inplace=True)
    hyper_dfs.append(hyper_df)

    size.viz.diagnostic_size_viz(
        samples, d, f'../../figures/size_diagnostics/')

# Concatenate and save the hyperparameters dataframes
hyper_df = pd.concat(hyper_dfs, sort=False)
summary_df = pd.concat(summary_dfs, sort=False)

# %%
hyper_df.to_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_samples.csv', index=False)
summary_df.to_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv', index=False)

# %%
