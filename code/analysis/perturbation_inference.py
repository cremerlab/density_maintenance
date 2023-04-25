# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
import size.viz
import seaborn as sns
import tqdm
import scipy.stats
cor, pal = size.viz.matplotlib_style()

# Compile the statistical model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/perturbation_inference_model.stan')

# %%
# ##############################################################################
# DATASET LOADING
# ##############################################################################
cal_data = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
brad_data = pd.read_csv(
    '../../data/protein_quantification/bradford_periplasmic_protein.csv')
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
biomass_data = pd.read_csv(
    '../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')
flow_data = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
growth_data = pd.read_csv(
    '../../data/summaries/summarized_growth_measurements.csv')
tot_prot = pd.read_csv(
    '../../data/literature/collated_total_protein_per_od.csv')

# %%
# ##############################################################################
# DATASET FILTERING
# ##############################################################################
# Correct bradford data that is out of bounds.
brad_data['od_600nm_true'] = brad_data['od_600nm']
brad_data.loc[brad_data['od_600nm'] >= 0.45, 'od_600nm_true'] = np.exp(
    1.26 * np.log(brad_data[brad_data['od_600nm'] >= 0.45]['od_600nm'].values) + 0.25)
brad_data = brad_data[~((brad_data['overexpression'] != 'none') & (
    brad_data['inducer_conc_ng_mL'] == 0))]

# Keep only the bradford data with more than two replicates
brad_data = pd.concat([d for _, d in brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']) if len(d) > 2], sort=False)
brad_data = brad_data[brad_data['strain'].isin(
    ['wildtype',  'malE-rbsB-fliC-KO', 'lpp14'])]
brad_data = brad_data[brad_data['overexpression'].isin(
    ['none', 'malE', 'rbsB', 'lacZ'])]

# Restrict size data
size_data = size_data[size_data['temperature_C'] == 37]
size_data = size_data[size_data['strain'].isin(
    ['wildtype', 'malE-rbsB-fliC-KO', 'lpp14'])]
size_data = pd.concat([d for _, d in size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']) if len(d) > 2], sort=False)
size_data = size_data[~((size_data['overexpression'] != 'none') & (
    size_data['inducer_conc'] == 0))]
size_data = size_data[size_data['overexpression'].isin(
    ['none', 'malE', 'rbsB', 'lacZ'])]
# Restrict flow data to only wildtype
flow_data = flow_data[flow_data['strain'] == 'wildtype']
flow_data = flow_data.groupby(
    ['date', 'carbon_source', 'run_no']).mean().reset_index()

# # Restrict growth data
growth_data = growth_data[growth_data['overexpression'].isin(['none', 'malE', 'rbsB', 'lacZ']) &
                          growth_data['strain'].isin(['wildtype', 'malE-rbsB-fliC-KO', 'lpp14'])]


# %%
# ##############################################################################
# DATA LABELING
# ##############################################################################
## SIZE INDEXING ###############################################################
# Add indexing to the size data
size_data['size_cond_idx'] = size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']).ngroup() + 1

## BRADFORD INDEXING ###########################################################
# Map size identifiers to the bradford conditions
brad_data['brad_mapper'] = 0
for g, d in size_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc', 'size_cond_idx']):
    brad_data.loc[(brad_data['strain'] == g[0]) &
                  (brad_data['carbon_source'] == g[1]) &
                  (brad_data['overexpression'] == g[2]) &
                  (brad_data['inducer_conc_ng_mL'] == g[3]),
                  'brad_mapper'] = g[-1]
brad_data = brad_data[brad_data['brad_mapper'] > 0]

# Filter, label, and transform bradford data
brad_data['cond_idx'] = brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']).ngroup() + 1
brad_data['conv_factor'] = brad_data['dilution_factor'] * \
    brad_data['extraction_volume_mL'] / \
    (brad_data['culture_volume_mL'])
brad_data['od_per_biomass'] = brad_data['od_595nm']


## GROWTH INDEXING #############################################################
growth_data['cond_idx'] = growth_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']).ngroup() + 1
for g, d in size_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc', 'size_cond_idx']):
    growth_data.loc[(growth_data['strain'] == g[0]) & (growth_data['carbon_source'] == g[1])
                    & (growth_data['overexpression'] == g[2]) & (growth_data['inducer_conc'] == g[3]),
                    'size_idx'] = g[-1]
growth_data.dropna(inplace=True)

## FLOW INDEXING ###############################################################
# Map size identifiers to the flow growth conditions.
flow_data['flow_mapper'] = 0
for g, d in size_data[(size_data['strain'] == 'wildtype') &
                      (size_data['overexpression'] == 'none') &
                      (size_data['inducer_conc'] == 0)
                      ].groupby(['carbon_source', 'size_cond_idx']):
    flow_data.loc[flow_data['carbon_source'] == g[0], 'flow_mapper'] = g[1]

# %%
data_dict = {

    'N_cal': len(cal_data),
    'concentration': cal_data['protein_conc_ug_ml'].values.astype(float),
    'cal_od': cal_data['od_595nm'].values.astype(float),
    'prot_frac': 0.55,

    'J_growth_cond': growth_data['cond_idx'].max(),
    'N_growth': len(growth_data),
    'growth_rates': growth_data['growth_rate_hr'].values.astype(float),
    'growth_idx': growth_data['cond_idx'].values.astype(int),
    'size_growth_idx': growth_data['size_idx'].values.astype(int),

    'N_size': len(size_data),
    'J_size_cond': size_data['size_cond_idx'].max(),
    'size_idx': size_data['size_cond_idx'].values.astype(int),
    'width': size_data['width_median'].values.astype(float),
    'length': size_data['length'].values.astype(float),
    'volume': size_data['volume'].values.astype(float),
    'peri_volume': size_data['periplasm_volume'].values.astype(float),
    'surface_area': size_data['surface_area'].values.astype(float),
    'surface_area_volume': size_data['surface_to_volume'].values.astype(float),
    'aspect_ratio': size_data['aspect_ratio'].values.astype(float),

    'N_brad': len(brad_data),
    'J_brad_cond': brad_data['cond_idx'].max(),
    'brad_idx': brad_data['cond_idx'].values.astype(int),
    'brad_mapper': brad_data['brad_mapper'].unique(),
    'brad_od595': brad_data['od_595nm'].values.astype(float),
    'brad_od600': brad_data['od_600nm'].values.astype(float),
    'conv_factor': brad_data['conv_factor'].values.astype(float),

    'N_biomass': len(biomass_data),
    'biomass': biomass_data['dry_mass_ug'],
    'N_prot': len(tot_prot),
    'total_protein': tot_prot['total_protein'],
    'total_protein_growth_rates': tot_prot['growth_rate_hr'],

    'N_flow': len(flow_data),
    'flow_mapper': flow_data['flow_mapper'].values.astype(int),
    'flow_events': flow_data['cells_per_biomass'].values.astype(float)
}


# %%
# Sample the posterior
_samples = model.sample(data_dict, adapt_delta=0.99)
samples = az.from_cmdstanpy(_samples)


# %%
size_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer_conc', 'size_cond_idx']
# Perform a KDE over the posteriors
shape_parameters = ['width_mu', 'length_mu', 'volume_mu', 'peri_volume_mu',
                    'alpha', 'growth_rate_rep', 'width_rep', 'length_rep',
                    'volume_rep', 'alpha_rep']
shape_post_kde = pd.DataFrame([])
for s in tqdm.tqdm(shape_parameters):
    post = samples.posterior[s].to_dataframe().reset_index()
    for g, d in size_data.groupby(size_groupby):
        for i, _g in enumerate(size_groupby[:-1]):
            post.loc[post[f'{s}_dim_0'] == g[-1]-1, _g] = g[i]
    minval = 0.01 * post[s].median()
    maxval = 3 * post[s].median()
    score_range = np.linspace(minval, maxval, 500)
    for g, d in post.groupby(size_groupby[:-1]):

        # Evaluate the kernel density
        kernel = scipy.stats.gaussian_kde(
            d[f'{s}'].values, bw_method='scott')
        kde = kernel(score_range)
        kde *= kde.sum()**-1

        # Set up the dataframe and store
        _df = pd.DataFrame(
            np.array([score_range, kde]).T, columns=['value', 'kde'])
        _df['parameter'] = s
        for i, _g in enumerate(size_groupby[:-1]):
            _df[_g] = g[i]
        shape_post_kde = pd.concat([shape_post_kde, _df], sort=False)
shape_post_kde.to_csv(
    '../../data/mcmc/perturbation_shape_posterior_kde.csv', index=False)

# %%
# Perform KDE over posterior of model params
brad_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer_conc_ng_mL', 'cond_idx']
model_params = ['phi_peri', 'm_peri', 'rho_peri', 'alpha']
model_post_kde = pd.DataFrame([])
for s in tqdm.tqdm(model_params):
    post = samples.posterior[s].to_dataframe().reset_index()
    for g, d in brad_data.groupby(brad_groupby):
        for i, _g in enumerate(brad_groupby[:-1]):
            post.loc[post[f'{s}_dim_0'] == g[-1]-1, _g] = g[i]
    minval = 0.01 * post[s].median()
    maxval = 3 * post[s].median()
    score_range = np.linspace(minval, maxval, 300)
    for g, d in post.groupby(brad_groupby[:-1]):

        # Evaluate the kernel density
        kernel = scipy.stats.gaussian_kde(
            d[f'{s}'].values, bw_method='scott')
        kde = kernel(score_range)
        kde *= kde.sum()**-1

        # Set up the dataframe and store
        _df = pd.DataFrame(
            np.array([score_range, kde]).T, columns=['value', 'kde'])
        _df['parameter'] = s
        for i, _g in enumerate(brad_groupby[:-1]):
            _df[_g] = g[i]
        model_post_kde = pd.concat([model_post_kde, _df], sort=False)
model_post_kde.to_csv(
    '../../data/mcmc/perturbation_parameter_posterior_kde.csv', index=False)


# %%
# Define the percentiles to keep
lower = [5, 12.5, 37.5, 50]
upper = [95, 87.5, 62.5, 50]
int_labels = ['95%', '75%', '25%', 'median']

# Process the parameters for the size inference
param_percs = pd.DataFrame([])
for i, p in enumerate(shape_parameters):
    p_df = samples.posterior[f'{p}'].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(p_df, p, f'{p}_dim_0',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=int_labels)
    for g, d in size_data.groupby(['strain', 'carbon_source', 'overexpression',
                                   'inducer_conc', 'size_cond_idx']):
        percs.loc[percs[f'{p}_dim_0'] == g[-1]-1, ['strain', 'carbon_source',
                                                   'overexpression', 'inducer_conc']] = g[:-1]
    percs.drop(columns=f'{p}_dim_0', inplace=True)
    param_percs = pd.concat([param_percs, percs], sort=False)

# Process the parameters for the protein quantities
params = ['phi_peri', 'm_peri', 'rho_peri']

for i, p in enumerate(params):
    p_df = samples.posterior[p].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(p_df, p, f'{p}_dim_0',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=int_labels)
    for g, d in brad_data.groupby(['strain', 'carbon_source', 'overexpression',
                                   'inducer_conc_ng_mL', 'cond_idx']):
        percs.loc[percs[f'{p}_dim_0'] == g[-1]-1, ['strain', 'carbon_source', 'overexpression',
                                                   'inducer_conc']] = g[:-1]
    percs.drop(columns=f'{p}_dim_0', inplace=True)
    param_percs = pd.concat([param_percs, percs])
param_percs.to_csv(
    '../../data/mcmc/perturbation_parameter_percentiles.csv', index=False)


# Process the parameters for the protein quantities
lam_df = samples.posterior['growth_mu'].to_dataframe().reset_index()
percs = size.viz.compute_percentiles(lam_df, 'growth_mu', 'growth_mu_dim_0',
                                     lower_bounds=lower,
                                     upper_bounds=upper,
                                     interval_labels=int_labels)
for g, d in growth_data.groupby(['strain', 'carbon_source', 'overexpression',
                                 'inducer_conc', 'cond_idx']):
    percs.loc[percs[f'growth_mu_dim_0'] == g[-1]-1, ['strain', 'carbon_source', 'overexpression',
                                                     'inducer_conc']] = g[:-1]
percs.drop(columns=f'growth_mu_dim_0', inplace=True)
percs.to_csv('../../data/mcmc/growth_parameter_percentiles.csv', index=False)
