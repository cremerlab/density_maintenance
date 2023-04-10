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
    '../../data/protein_quantification/bradford_periplasmic_protein.csv')
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
tot_prot = pd.read_csv('../../data/literature/collated_total_protein.csv')

# %%

# %%
# ##############################################################################
# DATASET FILTERING
# ##############################################################################
# Filter and aggregate mass spec
mass_spec_data = mass_spec_data[(mass_spec_data['periplasm'] == True)
                                ].groupby(['dataset_name',
                                           'growth_rate_hr',
                                           'condition']).sum().reset_index()

# Correct bradford data that is out of bounds.
brad_data['od_600nm_true'] = brad_data['od_600nm']
brad_data.loc[brad_data['od_600nm'] >= 0.45, 'od_600nm_true'] = np.exp(
    1.26 * np.log(brad_data[brad_data['od_600nm'] >= 0.45]['od_600nm'].values) + 0.25)
brad_data = brad_data[~((brad_data['overexpression'] != 'none') & (
    brad_data['inducer_conc_ng_mL'] == 0))]
# brad_data = brad_data[brad_data['od_600nm'] < 0.5]

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
J_growth_wt_idx = []
J_growth_size_idx = []
for g, d in size_data[(size_data['strain'] == 'wildtype') &
                      (size_data['overexpression'] == 'none') &
                      (size_data['inducer_conc'] == 0)].groupby(['carbon_source', 'size_cond_idx']):
    _growth = growth_data[(growth_data['strain'] == 'wildtype') &
                          (growth_data['carbon_source'] == g[0]) &
                          (growth_data['overexpression'] == 'none') &
                          (growth_data['inducer_conc'] == 0)]['cond_idx'].values[0]
    J_growth_wt_idx.append(_growth)
    J_growth_size_idx.append(g[-1])

# %%
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
    'growth_rates': growth_data['growth_rate_hr'].values.astype(float),
    'growth_cond_idx': growth_data['cond_idx'].values.astype(int),
    'J_growth_wt_idx': np.array(J_growth_wt_idx).astype(int),

    'size_growth_idx': growth_data['size_idx'].values.astype(int),
    'growth_size_idx': growth_data['size_idx'].values.astype(int),

    'J_size_wt': len(J_growth_size_idx),
    'J_size_wt_idx': np.array(J_growth_size_idx).astype(int),
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
    'N_prot': len(tot_prot),
    'total_protein': tot_prot['total_protein_ug_od600'],
    'total_protein_growth_rates': tot_prot['growth_rate_hr'],

    'N_mass_spec': len(mass_spec_data),
    'mass_fraction': mass_spec_data['mass_frac'],
    'mass_spec_growth_rate': mass_spec_data['growth_rate_hr'],

    'N_flow': len(flow_data),
    'flow_mapper': flow_data['flow_mapper'].values.astype(int),
    'cells_per_biomass': flow_data['cells_per_biomass'].values.astype(float)
}


# %%
# Sample the posterior
_samples = model.sample(data_dict, adapt_delta=0.99)
# adapt_delta=0.999, iter_sampling=2000)
samples = az.from_cmdstanpy(_samples)

# %%
size_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer_conc', 'size_cond_idx']
# Perform a KDE over the posteriors
shape_parameters = ['width_mu', 'length_mu', 'volume_mu', 'peri_volume_mu',
                    'surface_area_mu', 'surface_area_vol_mu', 'aspect_ratio_mu']
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

# %%
# Perform a KDE over the posteriors
lit_parameters = ['width_min', 'width_slope', 'length_min',
                  'length_slope', 'sav_min', 'sav_slope',
                  'alpha_min', 'alpha_slope',
                  'total_protein_min', 'total_protein_slope']

lit_post_kde = pd.DataFrame([])
for s in tqdm.tqdm(lit_parameters):
    post = samples.posterior[s].to_dataframe().reset_index()
    minval = 0.01 * post[s].median()
    maxval = 3 * post[s].median()
    score_range = np.linspace(minval, maxval, 500)
    # Evaluate the kernel density
    kernel = scipy.stats.gaussian_kde(post[f'{s}'].values, bw_method='scott')
    kde = kernel(score_range)
    kde *= kde.sum()**-1

    # Set up the dataframe and store
    _df = pd.DataFrame(
        np.array([score_range, kde]).T, columns=['value', 'kde'])
    _df['parameter'] = s
    lit_post_kde = pd.concat([lit_post_kde, _df], sort=False)
lit_post_kde.to_csv('../../data/mcmc/literature_model_posterior_kdes.csv')

# %%
# Perform KDE over posterior of model params
brad_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer_conc_ng_mL', 'cond_idx']
model_params = ['phi_M', 'prot_per_biomass_mu',
                'rho_peri', 'peri_prot_per_cell']
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
params = ['phi_M', 'prot_per_biomass_mu', 'rho_peri', 'peri_prot_per_cell']

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
param_percs.to_csv('../../data/mcmc/parameter_percentiles.csv', index=False)


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


# %%
singular_params = ['alpha_min', 'alpha_slope', 'width_min', 'width_slope',
                   'length_min', 'length_slope', 'sav_min', 'sav_slope',
                   'total_protein_min', 'total_protein_slope']
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
                                              interval_labels=int_labels)
singular_percs.drop(columns='idx', inplace=True)
singular_percs.to_csv(
    '../../data/mcmc/singular_parameter_percentiles.csv', index=False)
# %%
# Percentiles for mass spec data
params = ['mass_spec_phi_M', 'mass_spec_sav', 'mass_spec_widths', 'mass_spec_rho_peri',
          'mass_spec_peri_prot_per_cell', 'mass_spec_sa', 'mass_spec_tot_prot',
          'mass_spec_tot_prot_per_cell']

mass_spec_data['idx'] = np.arange(len(mass_spec_data))
mass_spec_df = pd.DataFrame([])
for i, p in enumerate(params):
    post = samples.posterior[p].to_dataframe().reset_index()
    percs = size.viz.compute_percentiles(post, p, f'{p}_dim_0',
                                         lower_bounds=lower,
                                         upper_bounds=upper,
                                         interval_labels=int_labels)
    for j, (g, d) in enumerate(percs.groupby(f'{p}_dim_0')):
        d['dataset_name'] = mass_spec_data['dataset_name'].values[j]
        d['growth_rate_hr'] = mass_spec_data['growth_rate_hr'].values[j]
        d['condition'] = mass_spec_data['condition'].values[j]
        d.drop(columns=f'{p}_dim_0', inplace=True)
        mass_spec_df = pd.concat([mass_spec_df, d], sort=False)
mass_spec_df.to_csv('../../data/mcmc/mass_spec_percentiles.csv', index=False)

# %%
# Growth rate dependence stuff
pairs = [['width_min', 'width_slope'],
         ['alpha_min', 'alpha_slope'],
         ['length_min', 'length_slope'],
         ['sav_min', 'sav_slope']]


lam_range = np.linspace(0.1, 2.5, 300)
perc_df = pd.DataFrame([])
for i, ell in enumerate(lam_range):
    for j, p in enumerate(pairs):
        if 'sav_min' == p[0]:
            pref = -1
        else:
            pref = 1
        # if 'alpha_min' == p[0]:
            # pref = 0
        pred_post = samples.posterior[p].to_dataframe().reset_index()
        pred_width = pred_post[p[0]] + pref * pred_post[p[1]] * ell
        _df = pd.DataFrame(np.array([pred_width]).T,
                           columns=['value'])
        _df['growth_rate_hr'] = ell
        _df['quantity'] = '_'.join(p[0].split('_')[:-1])
        percs = size.viz.compute_percentiles(_df, 'value',
                                             ['growth_rate_hr', 'quantity'],
                                             lower_bounds=lower,
                                             upper_bounds=upper,
                                             interval_labels=int_labels)
        perc_df = pd.concat([perc_df, percs])
perc_df.to_csv('../../data/mcmc/growth_rate_linear_relations.csv', index=False)

# %%
# ##############################################################################
# POSTERIOR PREDICTIVE CHECKS
# ##############################################################################
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
plt.savefig(
    '../../figures/mcmc/protein_diagnostics/one-shot_calibration_curve_ppc.pdf')
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
plt.savefig('../../figures/mcmc/protein_diagnostics/one-shot_bradford_ppc.pdf')

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
plt.savefig('../../figures/mcmc/protein_diagnostics/one-shot_biomass_ppc.pdf')
# %%
total_protein_ppc = samples.posterior.total_protein_rep.to_dataframe().reset_index()
for i in range(len(tot_prot)):
    total_protein_ppc.loc[total_protein_ppc['total_protein_rep_dim_0']
                          == i, 'growth_rate_hr'] = tot_prot['growth_rate_hr'].values[i]
    total_protein_ppc.loc[total_protein_ppc['total_protein_rep_dim_0']
                          == i, 'source'] = tot_prot['source'].values[i]
total_protein_perc = size.viz.compute_percentiles(total_protein_ppc, 'total_protein_rep',
                                                  ['growth_rate_hr', 'source'])

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
for i, (g, d) in enumerate(total_protein_perc.groupby(['interval'], sort=False)):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                    color=ppc_cmap[g], label='__nolegend__')

for g, d in tot_prot.groupby(['source']):
    ax.plot(d['growth_rate_hr'], d['total_protein_ug_od600'], 'o', label=g,
            ms=4)

ax.legend()
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('total protein [µg / OD$_{600nm}$ mL]')
plt.savefig(
    '../../figures/mcmc/protein_diagnostics/one-shot_total_protein_trend.pdf')

# %%
size_lam_percs = pd.DataFrame([])
params = ['widths_lit_rep', 'lengths_lit_rep', 'sav_lit_rep', 'alpha_lit_rep']
for i, p in enumerate(params):
    ppc = samples.posterior[p].to_dataframe().reset_index()
    for j in range(len(lit_size_data)):
        ppc.loc[ppc[f'{p}_dim_0'] == j,
                'growth_rate_hr'] = lit_size_data['growth_rate_hr'].values[j]
        ppc.loc[ppc[f'{p}_dim_0'] == j,
                'source'] = lit_size_data['source'].values[j]
    percs = size.viz.compute_percentiles(ppc, p, ['growth_rate_hr', 'source'])
    size_lam_percs = pd.concat([size_lam_percs, percs])

fig, ax = plt.subplots(4, 1, figsize=(4, 6), sharex=True)
axes = {'widths_lit_rep': ax[0], 'lengths_lit_rep': ax[1],
        'sav_lit_rep': ax[2], 'alpha_lit_rep': ax[3]}
for i, (g, d) in enumerate(size_lam_percs.groupby(['quantity', 'interval'], sort=False)):
    axes[g[0]].fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                            color=ppc_cmap[g[1]])

for g, d in lit_size_data.groupby(['source']):
    ax[0].plot(d['growth_rate_hr'], d['width_um'], 'o', label=g, ms=4)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], 'o', label=g, ms=4)
    ax[2].plot(d['growth_rate_hr'], d['surface_to_volume'], 'o', label=g, ms=4)
    ax[3].plot(d['growth_rate_hr'], d['length_um'] /
               d['width_um'], 'o', label=g, ms=4)
ax[0].set(title='literature width trend', xlabel='width [µm]')
ax[1].set(title='literature length trend', xlabel='length [µm]')
ax[2].set(title='literature SAV trend', xlabel='surface-to-volume [µm$^{-1}$]')
ax[3].set(title='literature aspect ratio', xlabel='aspect_ratio')
ax[0].set_ylim([0.5, 1.5])
ax[1].set_ylim([0, 7])
ax[2].set_ylim([2, 10])
ax[3].set_ylim([0, 7])
ax[3].set_xlabel('growth rate [hr$^{-1}$]')
plt.savefig('../../figures/mcmc/size_diagnostics/one-shot_literature_size_growth_trend.pdf',
            bbox_inches='tight')


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
plt.savefig('../../figures/mcmc/size_diagnostics/one-shot_size_ppc.pdf')

# %%
