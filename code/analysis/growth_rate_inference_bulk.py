# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import size.viz
import arviz
import tqdm
cor, pal = size.viz.matplotlib_style()

# Load the processed data
data = pd.read_csv(
    '../../data/growth_curves/wt_growth_measurements_processed.csv')

# Load and compile the inference model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/hierarchical_growth_inference.stan')

# %%
lam_df = pd.DataFrame([])
for g, d in data.groupby(['carbon_source', 'replicate']):
    # Assign grouping for hierarchical analysis
    d = d.copy()
    d['idx'] = (d.groupby(
        'date').ngroup() + 1).astype(int)

    # Define the data dictionary
    data_dict = {'J': d['idx'].max(),
                 'N': len(d),
                 'idx': d['idx'].values,
                 'elapsed_time': d['elapsed_time_hr'].values.astype(float),
                 'optical_density': d['od_600nm'].values.astype(float)}
    if data_dict['J'] == 1:
        print(f'Skipping sample {g} as there is only one replicate.')
        continue
    samples = model.sample(data=data_dict, iter_warmup=5000, iter_sampling=3000,
                           adapt_delta=0.99, show_progress=True)
    samples = arviz.from_cmdstanpy(samples)
    samples.to_netcdf(
        f'../../data/mcmc/wildtype_{g}_growth_inference_object.netcdf')

    # Make the diagnostic plots
    strain = d['strain'].values[0]
    carbon = g[0]
    temp = '37'
    oe = 'none'
    ind = 0
    metadata = {'strain': strain,
                'carbon': carbon,
                'temp': int(temp),
                'oe': oe,
                'ind': ind}
    size.viz.diagnostic_growth_viz(
        samples, d, metadata, '../../figures/mcmc/growth_diagnostics/')
    samples_df = samples.posterior['mu'].to_dataframe().reset_index()

    print('Processing growth rate....')
    lam = [d['mu'].values[0] for _, d in samples_df.groupby(['chain', 'draw'])]
    print('done!')
    _df = pd.DataFrame([])
    _df['replicate'] = g[1]
    _df['growth_rate_hr'] = lam
    _df['carbon_source'] = g[0]
    lam_df = pd.concat([lam_df, _df], sort=False)

# %%
# Compute the percentiles of the samples.
percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
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

# %%
