# %%
import numpy as np
import pandas as pd
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
    stan_file='./stan_models/one-shot_complete_inference.stan')

# %%
# ##############################################################################
# DATASET LOADING
# ##############################################################################
cal_data = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
brad_data = pd.read_csv(
    '../../data/protein_quantification/bradford_periplasmic_protein_v2.csv')
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
biomass_data = pd.read_csv(
    '../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')
flow_data = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
growth_data = pd.read_csv(
    '../../data/summaries/summarized_growth_measurements.csv')
lit_size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
mass_spec_data = pd.read_csv(
    '../../data/literature/compiled_mass_fractions.csv')

# ##############################################################################
# DATASET FILTERING
# ##############################################################################

# Filter and aggregate mass spec
mass_spec_data = mass_spec_data[(mass_spec_data['periplasm'] == True) &
                                (mass_spec_data['growth_rate_hr'] <= 0.8)
                                ].groupby(['dataset_name',
                                           'growth_rate_hr',
                                           'condition']).sum().reset_index()


# Correct bradford data that is out of bounds.
# TODO: (Maybe) Include into inference pipeline
brad_data['od_600nm_true'] = brad_data['od_600nm']
brad_data.loc[brad_data['od_600nm'] >= 0.45, 'od_600nm_true'] = np.exp(
    1.26 * np.log(brad_data[brad_data['od_600nm'] >= 0.45]['od_600nm'].values) + 0.25)
brad_data = brad_data[brad_data['od_600nm'] <= 0.5]

brad_data = brad_data[~((brad_data['overexpression'] != 'none') & (
    brad_data['inducer_conc_ng_mL'] == 0))]

# Keep only the bradford data with more than two replicates
brad_data = pd.concat([d for _, d in brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']) if len(d) > 2], sort=False)
brad_data = brad_data[brad_data['strain'].isin(
    ['wildtype',  'malE-rbsB-fliC-KO'])]

# Restrict size data
size_data = size_data[size_data['temperature_C'] == 37]
size_data = size_data[size_data['strain'].isin(
    ['wildtype', 'malE-rbsB-fliC-KO'])]
size_data = pd.concat([d for _, d in size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']) if len(d) > 2], sort=False)
size_data = size_data[~((size_data['overexpression'] != 'none') & (
    size_data['inducer_conc'] == 0))]

# Restrict flow data to only wildtype
flow_data = flow_data[flow_data['strain'] == 'wildtype']
flow_data = flow_data.groupby(
    ['date', 'carbon_source', 'run_no']).mean().reset_index()

# Restrict growth data
growth_data = growth_data[(growth_data['strain'] == 'wildtype') &
                          (growth_data['overexpression'] == 'none') &
                          (growth_data['inducer_conc'] == 0)]

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

# ##############################################################################
# DATA PREPARATION
# ##############################################################################
# Define the data dictionary
wt_brad = brad_data[(brad_data['strain'] == 'wildtype') &
                    (brad_data['overexpression'] == 'none') &
                    (brad_data['inducer_conc_ng_mL'] == 0)]

wt_size = size_data[(size_data['strain'] == 'wildtype') &
                    (size_data['overexpression'] == 'none') &
                    (size_data['inducer_conc'] == 0)]

data_dict = {

    'N_cal': len(cal_data),
    'concentration': cal_data['protein_conc_ug_ml'].values.astype(float),
    'cal_od': cal_data['od_595nm'].values.astype(float),
    'delta': size_data['delta'].unique()[0],
    'prot_frac': 0.55,

    'J_growth': growth_data['cond_idx'].max(),
    'N_growth': len(growth_data),
    'N_growth_lit': len(lit_size_data),
    'growth_rates_lit': lit_size_data['growth_rate_hr'].values.astype(float),
    'widths_lit': lit_size_data['width_um'].values.astype(float),
    'lengths_lit': lit_size_data['length_um'].values.astype(float),
    'growth_rates': growth_data['growth_rate_hr'].values.astype(float),
    'growth_cond_idx': growth_data['cond_idx'].values.astype(int),
    'growth_size_idx': growth_data['size_idx'].values.astype(int),

    'N_size': len(size_data),
    'N_size_wt': len(wt_size),
    'size_wt_idx': wt_size['size_cond_idx'].values.astype(int),
    'J_size_cond': size_data['size_cond_idx'].max(),
    'size_cond_idx': size_data['size_cond_idx'].values.astype(int),
    'width': size_data['width_median'].values.astype(float),
    'length': size_data['length'].values.astype(float),
    'volume': size_data['volume'].values.astype(float),
    'peri_volume': size_data['periplasm_volume'].values.astype(float),
    'surface_area': size_data['surface_area'].values.astype(float),
    'surface_area_volume': size_data['surface_to_volume'].values.astype(float),
    'aspect_ratio': size_data['aspect_ratio'].values.astype(float),

    'N_brad': len(brad_data),
    'N_brad_wt': len(wt_brad),
    'brad_wt_idx': wt_brad['cond_idx'].values.astype(int),
    'J_brad_cond': brad_data['cond_idx'].max(),
    'brad_cond_idx': brad_data['cond_idx'].values.astype(int),
    'brad_cond_mapper': brad_data['brad_mapper'].unique(),
    'brad_od595': brad_data['od_595nm'].values.astype(float),
    'brad_od600': brad_data['od_600nm'].values.astype(float),
    'conv_factor': brad_data['conv_factor'].values.astype(float),

    'N_biomass': len(biomass_data),
    'biomass': biomass_data['dry_mass_ug'],

    'N_mass_spec': len(mass_spec_data),
    'mass_fraction': mass_spec_data['mass_frac'],
    'mass_spec_growth_rate': mass_spec_data['growth_rate_hr'],

    'N_flow': len(flow_data),
    'flow_mapper': flow_data['flow_mapper'].values.astype(int),
    'cells_per_biomass': flow_data['cells_per_biomass'].values.astype(float)
}


# %%
# Sample the posterior
_samples = model.sample(data_dict,
                        adapt_delta=0.99, iter_sampling=2000)
samples = az.from_cmdstanpy(_samples)

# %%
size_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer_conc', 'size_cond_idx']
# Perform a KDE over the posteriors
shape_parameters = ['width_mu', 'length_mu', 'volume_mu', 'peri_volume_mu',
                    'surface_area_mu', 'surface_area_vol_mu', 'aspect_ratio_mu', ]
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
shape_post_kde.to_csv('../../data/mcmc/shape_posterior_kde.csv', index=False)

# Perform KDE over posterior of model params
brad_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer_conc_ng_mL', 'cond_idx']
model_params = ['phi_M', 'prot_per_biomass_mu', 'prot_per_cell']
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
model_post_kde.to_csv('../../data/mcmc/model_posterior_kde.csv', index=False)

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
model_post_kde.to_csv('../../data/mcmc/model_posterior_kde.csv', index=False)


# %%
# Define the percentiles to keep
lower = [5, 12.5, 37.5, 50]
upper = [95, 87.5, 62.5, 50]
labels = ['95%', '75%', '25%', 'median']
# Process the parameters for the size inference
param_percs = pd.DataFrame([])
for i, p in enumerate(shape_parameters):
    p_df = samples.posterior[f'{p}'].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(p_df, p, f'{p}_dim_0',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=labels)
    for g, d in size_data.groupby(['strain', 'carbon_source', 'overexpression',
                                   'inducer_conc', 'size_cond_idx']):
        percs.loc[percs[f'{p}_dim_0'] == g[-1]-1, ['strain', 'carbon_source',
                                                   'overexpression', 'inducer_conc']] = g[:-1]
    percs.drop(columns=f'{p}_dim_0', inplace=True)
    param_percs = pd.concat([param_percs, percs], sort=False)

# Process the parameters for the protein quantities
params = ['phi_M', 'prot_per_biomass_mu', 'N_cells', 'prot_per_cell']
for i, p in enumerate(params):
    p_df = samples.posterior[p].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(p_df, p, f'{p}_dim_0',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=labels)
    for g, d in brad_data.groupby(['strain', 'carbon_source', 'overexpression',
                                   'inducer_conc_ng_mL', 'cond_idx']):
        percs.loc[percs[f'{p}_dim_0'] == g[-1]-1, ['strain', 'carbon_source', 'overexpression',
                                                   'inducer_conc']] = g[:-1]
    percs.drop(columns=f'{p}_dim_0', inplace=True)
    param_percs = pd.concat([param_percs, percs])
param_percs.to_csv('../../data/mcmc/parameter_percentiles.csv', index=False)

# %%
singular_params = ['Lambda', 'alpha', 'phi_max', 'width_min',
                   'avg_mass_spec_k', 'slope', 'intercept']
singular_percs = samples.posterior[singular_params].to_dataframe(
).reset_index()
singular_kde = pd.DataFrame([])
for p in singular_params:
    kernel = scipy.stats.gaussian_kde(singular_percs[p], bw_method='scott')
    score_range = np.linspace(
        0.01 * singular_percs[p].median(), 3 * singular_percs[p].median(), 300)
    kde = kernel(score_range)
    kde *= kde.sum()**-1
    _df = pd.DataFrame(
        np.array([score_range, kde]).T, columns=['value', 'kde'])
    _df['parameter'] = p
    singular_kde = pd.concat([singular_kde, _df])
singular_kde.to_csv('../../data/mcmc/singular_posterior_kde.csv', index=False)
singular_percs['idx'] = 0
singular_percs = size.viz.compute_percentiles(singular_percs, singular_params, 'idx',
                                              lower_bounds=lower,
                                              upper_bounds=upper,
                                              interval_labels=labels)
singular_percs.drop(columns='idx', inplace=True)
singular_percs.to_csv(
    '../../data/mcmc/singular_parameter_percentiles.csv', index=False)

# %%
pred_post = samples.posterior[['Lambda', 'phi_max', 'width_min', 'intercept', 'slope',
                               'width_slope', 'alpha']].to_dataframe().reset_index()


lam_range = np.linspace(0.1, 2.5, 300)
perc_df = pd.DataFrame([])
for i, ell in enumerate(lam_range):
    pred_width = pred_post['width_min'] + pred_post['width_slope'] * ell
    pred_sav = 12 * pred_post['alpha'] / \
        (pred_width * (3 * pred_post['alpha'] - 1))
    _df = pd.DataFrame(np.array([pred_width, pred_sav]).T,
                       columns=['width_um', 'sav_inv_um'])
    _df['growth_rate_hr'] = ell
    percs = size.viz.compute_percentiles(_df, ['width_um', 'sav_inv_um'],
                                         'growth_rate_hr',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=labels)
    perc_df = pd.concat([perc_df, percs])
perc_df.to_csv('../../data/mcmc/dimensions_v_growth_rate.csv', index=False)

# %%
mass_spec_conv = pd.DataFrame([])
for val in ['mass_spec_widths', 'mass_spec_sav']:
    _post = samples.posterior[val].to_dataframe().reset_index()
    _post = _post.groupby([f'{val}_dim_0']).median().reset_index()
    _post['dataset_name'] = mass_spec_data['dataset_name']
    _post['growth_rate_hr'] = mass_spec_data['growth_rate_hr']
    _post['phi_M'] = data_dict['mass_fraction'] * data_dict['prot_frac']
    _post = _post[['dataset_name', 'growth_rate_hr', 'phi_M', f'{val}']]
    _post.rename(columns={f'{val}': 'value'}, inplace=True)
    if 'widths' in val:
        val = val[:-1]
    _post['quantity'] = val.split('_')[-1]

    mass_spec_conv = pd.concat([mass_spec_conv, _post])
mass_spec_conv.to_csv(
    '../../data/mcmc/converted_mass_spec_medians.csv', index=False)
# %%
# Compute the two predictions and the percentiles
# phi_range = np.linspace(0.001, 0.1, 200)
# perc_df = pd.DataFrame([])
# for i, p in enumerate(tqdm.tqdm(phi_range)):
#     taylor = pred_post['width_min'] + pred_post['Lambdak'] * ((p - pred_post['phi_star'])/pred_post['phi_star'] ** 2) *\
#         ((p - pred_post['phi_star'])/pred_post['phi_star'] - 1)
#     simple = pred_post['Lambda'] * \
#         (pred_post['avg_rho_ratio'] * (1 / p - 1) + 1)

#     _df = pd.DataFrame(np.array([taylor, simple]).T,
#                        columns=['width_taylor', 'width_simple'])
#     _df['phi_M'] = p
#     percs = size.viz.compute_percentiles(_df, ['width_taylor', 'width_simple'], 'phi_M',
#                                          lower_bounds=lower,
#                                          upper_bounds=upper,
#                                          interval_labels=labels)
#     perc_df = pd.concat([perc_df, percs], sort=False)
# perc_df.to_csv('../../data/mcmc/predicted_scaling_lppwt.csv', index=False)

# %
# Compute the two predictions and the percentiles
w_range = np.linspace(0.45, 1, 200)
perc_df = pd.DataFrame([])
for i, w in enumerate(tqdm.tqdm(w_range)):
    pred = pred_post['intercept'] + pred_post['slope'] * w

    _df = pd.DataFrame(np.array([pred]).T,
                       columns=['phi_pred'])
    _df['width'] = w
    percs = size.viz.compute_percentiles(_df, 'phi_pred', 'width',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=labels)
    perc_df = pd.concat([perc_df, percs], sort=False)
perc_df.to_csv('../../data/mcmc/predicted_phi_scaling_lppwt.csv', index=False)

# %%
# Plot the ppc for the calibration data
cal_ppc = samples.posterior.od595_calib_rep.to_dataframe().reset_index()
for conc, n in zip(cal_data['protein_conc_ug_ml'], np.arange(len(cal_data))):
    cal_ppc.loc[cal_ppc['od595_calib_rep_dim_0'] == n, 'conc'] = conc
cal_percs = size.viz.compute_percentiles(cal_ppc, 'od595_calib_rep', 'conc')

fig, ax = plt.subplots(1, 1)
cmap = sns.color_palette('Greys', len(cal_percs['interval'].unique()) + 1)
ppc_cmap = {i: c for i, c in zip(cal_percs['interval'].unique(), cmap)}
for i, (g, d) in enumerate(cal_percs.groupby(['interval'], sort=False)):
    ax.fill_between(d['conc'].values, d['lower'].values, d['upper'].values,
                    color=ppc_cmap[g], zorder=i+1)
ax.plot(cal_data['protein_conc_ug_ml'], cal_data['od_595nm'], 'o',
        color=cor['primary_red'], ms=4, zorder=i+1)
ax.set_xlabel('protein standard concentration [µg / mL]')
ax.set_ylabel('OD$_{595nm}$')

# %%
# Compute the percentiles
prot_ppc_df = samples.posterior.od595_brad_rep.to_dataframe().reset_index()

for dim, idx in zip(np.arange(len(brad_data)), brad_data['cond_idx'].values):
    prot_ppc_df.loc[prot_ppc_df['od595_brad_rep_dim_0']
                    == dim, 'cond_idx'] = idx

prot_ppc_percs = size.viz.compute_percentiles(
    prot_ppc_df, 'od595_brad_rep', 'cond_idx')

fig, ax = plt.subplots(1, 1, figsize=(4, 6))

for i, (g, d) in enumerate(prot_ppc_percs.groupby(['cond_idx', 'interval'], sort=False)):
    ax.hlines(g[0], d['lower'], d['upper'], lw=5,
              color=ppc_cmap[g[1]], zorder=i+1)


for g, d in brad_data.groupby(['cond_idx']):
    ax.plot(d['od_per_biomass'], np.ones(len(d)) * g + np.random.normal(0, 0.05, len(d)),
            'o', color=cor['primary_red'], ms=4, zorder=i+1)

labels = []
for g, d in brad_data.groupby(['cond_idx', 'strain', 'carbon_source',
                               'overexpression', 'inducer_conc_ng_mL']):
    labels.append(g[1:])
ax.set_yticks(brad_data['cond_idx'].unique())
ax.set_yticklabels(labels)
ax.set_xlabel('OD$_{595nm}$')

# %%
# Biomass ppc
biomass_ppc = samples.posterior.biomass_rep.to_dataframe()
biomass_ppc['idx'] = 1
biomass_percs = size.viz.compute_percentiles(biomass_ppc, 'biomass_rep', 'idx')

fig, ax = plt.subplots(1, 1, figsize=(4, 1))

for i, (g, d) in enumerate(biomass_percs.groupby(['interval'], sort=False)):
    ax.hlines(1, d['lower'], d['upper'], lw=15, zorder=i+1, color=ppc_cmap[g])


ax.plot(biomass_data['dry_mass_ug'], np.ones(
    len(biomass_data)), 'o', color=cor['primary_green'], zorder=i+1)
ax.set_yticks([])
ax.set_xlabel('dry mass [µg / OD$_{600nm}$ mL]')


# %%
fig, ax = plt.subplots(3, 3, figsize=(6, 12), sharey=True)
props = ['width', 'length', 'volume', 'peri_volume',
         'surface_area', 'surface_area_vol', 'aspect_ratio']
for p, a in zip(props, ax.ravel()):
    a.set_xlabel(p)

loc = {k: ax.ravel()[i] for i, k in enumerate(props)}
for p in props:
    _ppc = samples.posterior[f'{p}_rep'].to_dataframe().reset_index()
    for dim, idx in zip(np.arange(len(size_data)), size_data['size_cond_idx'].values):
        _ppc.loc[_ppc[f'{p}_rep_dim_0'] == dim, 'idx'] = idx
    _perc = size.viz.compute_percentiles(_ppc, f'{p}_rep', 'idx')
    for i, (g, d) in enumerate(_perc.groupby(['idx', 'interval'], sort=False)):
        loc[p].hlines(g[0], d['lower'], d['upper'], lw=2, color=ppc_cmap[g[1]],
                      zorder=i+1)

for g, d in size_data.groupby(['size_cond_idx']):
    _ones = np.ones(len(d))
    ax[0, 0].plot(d['width_median'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)
    ax[0, 1].plot(d['length'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)
    ax[0, 2].plot(d['volume'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)
    ax[1, 0].plot(d['periplasm_volume'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)
    ax[1, 1].plot(d['surface_area'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)
    ax[1, 2].plot(d['surface_to_volume'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)
    ax[2, 0].plot(d['aspect_ratio'], _ones * g +
                  np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
                  ms=3, zorder=1000)

labels = [g[1:] for g, _ in size_data.groupby(
    ['size_cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])]
_ = ax[0, 0].set_yticks(size_data['size_cond_idx'].unique())
_ = ax[1, 0].set_yticks(size_data['size_cond_idx'].unique())
_ = ax[0, 0].set_yticklabels(labels)
_ = ax[1, 0].set_yticklabels(labels)


# # %%

# # Link the axes
# brad_idx = {g[1:]: g[0] for g, _ in brad_data.groupby(
#     ['cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL'])}
# size_loc = {g[0]: brad_idx[g[1:]] for g, _ in size_data.groupby(
#     ['size_cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc']) if g[1:] in brad_idx.keys()}
# width_ppc = samples.posterior.width_mu.to_dataframe().reset_index()

# frac_ppc = samples.posterior.phi_M.to_dataframe().reset_index()

# for g, _ in brad_data.groupby(['cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']):
#     frac_ppc.loc[frac_ppc['phi_M_dim_0'] == g[0]-1, [
#         'strain', 'carbon_source', 'overexpression', 'inducer_conc']] = g[1:]
#     frac_ppc.loc[frac_ppc['phi_M_dim_0']
#                  == g[0]-1, 'link_idx'] = brad_idx[g[1:]]
# frac_perc = size.viz.compute_percentiles(frac_ppc, 'phi_M', [
#                                          'link_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])

# for g, _ in size_data.groupby(['size_cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc']):
#     if g[0] in size_loc:
#         width_ppc.loc[width_ppc['width_mu_dim_0'] == g[0]-1, [
#             'strain', 'carbon_source', 'overexpression', 'inducer_conc']] = g[1:]
#         width_ppc.loc[width_ppc['width_mu_dim_0']
#                       == g[0]-1, 'link_idx'] = size_loc[g[0]]

# # width_ppc.dropna(inplace=True)
# width_perc = size.viz.compute_percentiles(width_ppc, 'width_mu', [
#                                           'link_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])


# # %%
# # Mass spec data
# mass_fracs = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
# mass_fracs = mass_fracs[mass_fracs['periplasm']]
# dry_frac = 0.3
# prot_frac = 0.55
# density = 1.1

# # %%
# # Do the proper classification
# # genes = pd.read_csv('../../data/literature/genes_classification_all.csv')
# # _genes = genes[genes['location'].isin(['IM', 'OM', 'PE', 'LPO'])]
# # mass_fracs = mass_fracs[mass_fracs['gene_name'].isin(_genes['gene'].unique())]
# growth = pd.read_csv('../../data/summaries/summarized_growth_measurements.csv')
# growth = growth[(growth['strain'] == 'wildtype') & (
#     growth['overexpression'] == 'none') & (growth['inducer_conc_ng_ml'] == 0)]
# growth = growth.groupby(['carbon_source']).mean().reset_index()
# wt_size = size_data[(size_data['strain'] == 'wildtype') & (
#     size_data['overexpression'] == 'none') & (size_data['inducer_conc'] == 0)]
# for g, d in growth.groupby(['carbon_source', 'growth_rate_hr']):
#     wt_size.loc[wt_size['carbon_source'] == g[0], 'growth_mu'] = g[1]

# # Determine the simple relations
# w_popt = scipy.stats.linregress(
#     wt_size['growth_mu'], wt_size['width_median'])
# ell_popt = scipy.stats.linregress(
#     wt_size['growth_mu'], np.log(wt_size['length']))
# peri_vol_popt = scipy.stats.linregress(
#     wt_size['growth_mu'], np.log(wt_size['periplasm_volume']))
# vol_popt = scipy.stats.linregress(
#     wt_size['growth_mu'], np.log(wt_size['volume']))
# sav_popt = scipy.stats.linregress(
#     wt_size['growth_mu'], np.log(wt_size['surface_to_volume']))

# # Compute the periplasmic protein density
# mass_fracs['width'] = w_popt[0] * \
#     mass_fracs['growth_rate_hr'] + w_popt[1]
# mass_fracs['length'] = np.exp(
#     ell_popt[0] * mass_fracs['growth_rate_hr'] + ell_popt[1])
# mass_fracs['peri_vol'] = np.exp(
#     peri_vol_popt[0] * mass_fracs['growth_rate_hr'] + peri_vol_popt[1])
# mass_fracs['volume'] = np.exp(
#     vol_popt[0] * mass_fracs['growth_rate_hr'] + vol_popt[1])
# mass_fracs['sav'] = np.exp(
#     sav_popt[0] * mass_fracs['growth_rate_hr'] + sav_popt[1])
# # mass_fracs['peri_volume'] = size.analytical.surface_area(mass_fracs['length'], mass_fracs['width']) * 0.025
# mass_fracs['tot_protein'] = density * \
#     dry_frac * prot_frac * mass_fracs['volume']
# mass_fracs['peri_protein'] = mass_fracs['mass_frac'] * \
#     mass_fracs['tot_protein']
# mass_fracs['rho_peri'] = (
#     mass_fracs['peri_protein'] * 1E3) / mass_fracs['peri_vol']
# mass_fracs['biomass_frac'] = mass_fracs['peri_protein'] / \
#     (density * dry_frac * mass_fracs['volume'])
# mass_fracs = mass_fracs.groupby(
#     ['dataset_name', 'condition', 'growth_rate_hr', 'width']).sum().reset_index()
