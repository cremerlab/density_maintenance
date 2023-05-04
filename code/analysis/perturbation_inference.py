# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner.corner as corner
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
brad_data = pd.read_csv(
    '../../data/processed/labeled_protein_measurements.csv')
growth_data = pd.read_csv(
    '../../data/processed/labeled_growth_measurements.csv')
size_data = pd.read_csv('../../data/processed/labeled_size_measurements.csv')
biuret_data = pd.read_csv(
    '../../data/processed/labeled_total_protein_measurements.csv')
flow_data = pd.read_csv('../../data/processed/labeled_flow_measurements.csv')

brad_cal = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
biuret_cal = pd.read_csv(
    '../../data/protein_quantification/biuret_calibration_curve.csv')

brad_data['conv_factor'] = brad_data['dilution_factor'] * \
    brad_data['extraction_volume_mL'] / \
    (brad_data['culture_volume_mL'])
# %%
# Assemble the data dictionary
data_dict = {
    'N_size': len(size_data),
    'N_brad': len(brad_data),
    'N_biuret': len(biuret_data),
    'N_flow': len(flow_data),
    'N_growth': len(growth_data),
    'N_biuret_cal': len(biuret_cal),
    'N_brad_cal': len(brad_cal),
    'N_biuret_cal': len(biuret_cal),

    'J_cond': size_data['perturbation_idx'].max(),
    'J_brad_cond': brad_data['perturbation_idx'].max(),
    'J_flow_cond': len(flow_data['carbon_idx'].unique()),

    'width': size_data['width_median'].values,
    'length': size_data['length'].values,
    'volume': size_data['volume'].values,
    'peri_volume': size_data['periplasm_volume'].values,
    'size_idx': size_data['perturbation_idx'].values,

    'growth_rates': growth_data['growth_rate_hr'].values,
    'growth_idx': growth_data['perturbation_idx'].values,

    'flow_events': flow_data['cells_per_biomass'].values,
    'flow_idx':  flow_data['cond_idx'].values,
    'flow_mapper': flow_data['wt_carbon_idx'].values,

    'brad_od595': brad_data['od_595nm'].values,
    'brad_od600': brad_data['od_600nm'].values,
    'brad_conv_factor': brad_data['conv_factor'].values,
    'brad_cal_conc': brad_cal['protein_conc_ug_ml'].values,
    'brad_cal_od': brad_cal['od_595nm'].values,
    'brad_idx': brad_data['perturbation_idx'].values,
    'brad_mapper': brad_data['perturbation_idx'].unique(),
    'brad_flow_mapper': brad_data.groupby(['perturbation_idx'])['carbon_idx'].min(),

    'biuret_od555': biuret_data['od_555nm'].values,
    'biuret_od600': biuret_data['od_600nm'].values,
    'biuret_conv_factor': biuret_data['conv_factor'].values,
    'biuret_cal_od': biuret_cal['od_555nm'].values,
    'biuret_cal_conc': biuret_cal['protein_conc_ug_ml'].values,
    'biuret_mapper': biuret_data['wt_carbon_idx'].values
}

# |, show_console=True)
# , show_console=True)
_samples = model.sample(data=data_dict, adapt_delta=0.95)
samples = az.from_cmdstanpy(_samples)

# %%
# ##############################################################################
# POSTERIOR PREDICTIVE CHECKS
# ##############################################################################
# Plot the ppc for the calibration data

##########
# BRADFORD
##########
cal_ppc = samples.posterior.od595_cal_rep.to_dataframe().reset_index()
for conc, n in zip(brad_cal['protein_conc_ug_ml'], np.arange(len(brad_cal))):
    cal_ppc.loc[cal_ppc['od595_cal_rep_dim_0'] == n, 'conc'] = conc
cal_percs = size.viz.compute_percentiles(cal_ppc, 'od595_cal_rep', 'conc')

fig, ax = plt.subplots(1, 1)
cmap = sns.color_palette('Greys', len(cal_percs['interval'].unique()) + 1)
ppc_cmap = {i: c for i, c in zip(cal_percs['interval'].unique(), cmap)}
for i, (g, d) in enumerate(cal_percs.groupby(['interval'], sort=False)):
    ax.fill_between(d['conc'].values, d['lower'].values, d['upper'].values,
                    color=ppc_cmap[g], zorder=i+1)
ax.plot(brad_cal['protein_conc_ug_ml'], brad_cal['od_595nm'], 'o',
        color=cor['primary_red'], ms=4, zorder=i+1)
ax.set_xlabel('protein standard concentration [µg / mL]')
ax.set_ylabel('OD$_{595nm}$')
plt.savefig(
    '../../figures/mcmc/protein_diagnostics/perturbation_inference_bradford_calibration_curve_ppc.pdf')

##########
# BIURET
##########
cal_ppc = samples.posterior.biuret_cal_rep.to_dataframe().reset_index()
for conc, n in zip(biuret_cal['protein_conc_ug_ml'], np.arange(len(brad_cal))):
    cal_ppc.loc[cal_ppc['biuret_cal_rep_dim_0'] == n, 'conc'] = conc
cal_percs = size.viz.compute_percentiles(cal_ppc, 'biuret_cal_rep', 'conc')

fig, ax = plt.subplots(1, 1)
cmap = sns.color_palette('Greys', len(cal_percs['interval'].unique()) + 1)
ppc_cmap = {i: c for i, c in zip(cal_percs['interval'].unique(), cmap)}
for i, (g, d) in enumerate(cal_percs.groupby(['interval'], sort=False)):
    ax.fill_between(d['conc'].values, d['lower'].values, d['upper'].values,
                    color=ppc_cmap[g], zorder=i+1)
ax.plot(biuret_cal['protein_conc_ug_ml'], biuret_cal['od_555nm'], 'o',
        color=cor['primary_red'], ms=4, zorder=i+1)
ax.set_xlabel('protein standard concentration [µg / mL]')
ax.set_ylabel('OD$_{555nm}$')
plt.savefig(
    '../../figures/mcmc/protein_diagnostics/perturbation_inference_biuret_calibration_curve_ppc.pdf')

# %%
# Compute the percentiles
prot_ppc_df = samples.posterior.od595_meas_rep.to_dataframe().reset_index()

for dim, idx in zip(np.arange(len(brad_data)), brad_data['perturbation_idx'].values):
    prot_ppc_df.loc[prot_ppc_df['od595_meas_rep_dim_0']
                    == dim, 'perturbation_idx'] = idx

prot_ppc_percs = size.viz.compute_percentiles(
    prot_ppc_df, 'od595_meas_rep', 'perturbation_idx')

fig, ax = plt.subplots(1, 1, figsize=(4, 6))

for i, (g, d) in enumerate(prot_ppc_percs.groupby(['perturbation_idx', 'interval'], sort=False)):
    ax.hlines(g[0], d['lower'], d['upper'], lw=5,
              color=ppc_cmap[g[1]], zorder=i+1)


for g, d in brad_data.groupby(['perturbation_idx']):
    ax.plot(d['od_595nm'], np.ones(len(d)) * g + np.random.normal(0, 0.05, len(d)),
            'o', color=cor['primary_red'], ms=4, zorder=i+1)

labels = []
for g, d in brad_data.groupby(['perturbation_idx', 'strain', 'carbon_source',
                               'overexpression', 'inducer_conc']):
    labels.append(g)
ax.set_yticks(np.arange(len(labels))+1)
ax.set_yticklabels(labels)
ax.set_xlabel('OD$_{595nm}$')
plt.savefig(
    '../../figures/mcmc/protein_diagnostics/perturbation_inference_bradford_ppc.pdf')

#  %%
total_protein_ppc = samples.posterior.total_protein_od555_rep.to_dataframe().reset_index()
for i in range(len(biuret_data)):
    total_protein_ppc.loc[total_protein_ppc['total_protein_od555_rep_dim_0']
                          == i, 'growth_rate_hr'] = biuret_data['growth_rate_hr'].values[i]
total_protein_perc = size.viz.compute_percentiles(total_protein_ppc, 'total_protein_od555_rep',
                                                  ['growth_rate_hr'])

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
for i, (g, d) in enumerate(total_protein_perc.groupby(['interval'], sort=False)):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'],
                    color=ppc_cmap[g], label='__nolegend__')

ax.plot(biuret_data['growth_rate_hr'], biuret_data['od_555nm'],
        'o', label=g, ms=4, color=cor['primary_red'])

ax.legend()
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('OD$_{555nm}$')
plt.savefig(
    '../../figures/mcmc/protein_diagnostics/perturbation_inference_total_protein_trend.pdf')
# %% Flow PPC
fig, ax = plt.subplots(1, 1, figsize=(4, 6))
ax.set_xscale('log')

ppc = samples.posterior.flow_rep.to_dataframe().reset_index()
for i in range(len(flow_data)):
    _d = flow_data.iloc[i]
    ppc.loc[ppc['flow_rep_dim_0'] == i, 'carbon_source'] = _d['carbon_source']
perc = size.viz.compute_percentiles(ppc, 'flow_rep', 'carbon_source')

labels = []
locs = {}
for i, (g, d) in enumerate(flow_data.groupby(['carbon_source'])):
    locs[g] = i
    ax.plot(d['cells_per_biomass'], locs[g] + np.random.normal(0, 0.05, len(d)), 'o',
            color=cor['primary_red'], zorder=1000)
    labels.append(g)


for i, (g, d) in enumerate(perc.groupby(['carbon_source', 'interval'], sort=False)):
    ax.hlines(locs[g[0]], d['lower'], d['upper'], lw=5, color=ppc_cmap[g[1]])

ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels)
ax.set_xlabel('cells per OD$_{600}$ mL')
# %%
fig, ax = plt.subplots(1, 5, figsize=(12, 4), sharey=True)
# ax = ax.ravel()
props = ['width', 'length', 'volume', 'peri_volume', 'alpha']
for p, a in zip(props, ax.ravel()):
    a.set_xlabel(p)

loc = {k: ax[i] for i, k in enumerate(props)}
for p in props:
    _ppc = samples.posterior[f'{p}_rep'].to_dataframe().reset_index()

    for idx in size_data['perturbation_idx'].unique():
        _ppc.loc[_ppc[f'{p}_rep_dim_0'] == idx - 1, 'idx'] = idx

    _perc = size.viz.compute_percentiles(_ppc, f'{p}_rep', 'idx')
    for i, (g, d) in enumerate(_perc.groupby(['idx', 'interval'], sort=False)):

        loc[p].hlines(g[0], d['lower'], d['upper'], lw=2, color=ppc_cmap[g[1]],
                      zorder=i+1)

labels = []
idx = 1
for g, d in size_data.groupby(['perturbation_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc']):
    labels.append(g[1:])
    _ones = np.ones(len(d))
    ax[0].plot(d['width_median'], _ones * idx +
               np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
               ms=3, zorder=1000)
    ax[1].plot(d['length'], _ones * idx +
               np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
               ms=3, zorder=1000)
    ax[2].plot(d['volume'], _ones * idx +
               np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
               ms=3, zorder=1000)
    ax[3].plot(d['periplasm_volume'], _ones * idx +
               np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
               ms=3, zorder=1000)
    # ax[1, 1].plot(d['surface_area'], _ones * g +
    #   np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
    #   ms=3, zorder=1000)
    # ax[1, 2].plot(d['surface_to_volume'], _ones * g +
    #   np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
    #   ms=3, zorder=1000)
    ax[4].plot(d['aspect_ratio'], _ones * idx +
               np.random.normal(0, 0.05, len(d)), 'o', color=cor['primary_red'],
               ms=3, zorder=1000)
    idx += 1

_ = ax[0].set_yticks(1 + np.arange(len(labels)))
_ = ax[3].set_yticks(1 + np.arange(len(labels)))
_ = ax[0].set_yticklabels(labels)
_ = ax[3].set_yticklabels(labels)
plt.savefig(
    '../../figures/mcmc/size_diagnostics/perturbation_inference_size_ppc.pdf')

# %%
# Visualize growth rate ppcs.
fig, ax = plt.subplots(1, 1, figsize=(4, 7))
inds = {g: i for i, (g, _) in enumerate(growth_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer', 'inducer_conc', 'temperature']))}
for g, d in growth_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer', 'inducer_conc', 'temperature']):
    ax.plot(d['growth_rate_hr'], inds[g] + np.random.normal(0, 0.1, len(d)), 'o', color=cor['primary_red'],
            ms=3)
ppcs = samples.posterior.growth_rate_rep.to_dataframe().reset_index()
for i in range(len(growth_data)):
    _d = growth_data.iloc[i]
    ppcs.loc[ppcs['growth_rate_rep_dim_0'] == i, 'strain'] = _d['strain']
    ppcs.loc[ppcs['growth_rate_rep_dim_0'] == i,
             'carbon_source'] = _d['carbon_source']
    ppcs.loc[ppcs['growth_rate_rep_dim_0'] == i,
             'overexpression'] = _d['overexpression']
    ppcs.loc[ppcs['growth_rate_rep_dim_0'] ==
             i, 'inducer_conc'] = _d['inducer_conc']
    ppcs.loc[ppcs['growth_rate_rep_dim_0'] ==
             i, 'inducer'] = _d['inducer']
    ppcs.loc[ppcs['growth_rate_rep_dim_0'] ==
             i, 'temperature'] = _d['temperature']

percs = size.viz.compute_percentiles(ppcs, 'growth_rate_rep', [
                                     'strain', 'carbon_source', 'overexpression', 'inducer', 'inducer_conc', 'temperature'])

for g, d in percs.groupby(['strain', 'carbon_source', 'overexpression', 'inducer', 'inducer_conc', 'temperature', 'interval'], sort=False):
    ax.hlines(inds[g[:-1]], d['lower'], d['upper'],
              lw=10, color=ppc_cmap[g[-1]], zorder=1)

ax.set_yticks(list(inds.values()))
ax.set_yticklabels(list(inds.keys()))
ax.set_xlabel('growth rate [hr$^{-1}$]')
# %%
# ##############################################################################
# PARAMETER SUMMARIZATION
# ##############################################################################
size_groupby = ['strain', 'carbon_source',
                'overexpression', 'inducer', 'inducer_conc', 'temperature', 'perturbation_idx']
# Perform a KDE over the posteriors
shape_parameters = ['width_mu', 'length_mu', 'volume_mu', 'peri_volume_mu',
                    'alpha_mu', 'alpha_rep', 'growth_rate_rep', 'width_rep', 'length_rep',
                    'volume_rep', 'peri_volume_rep']
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
                'overexpression', 'inducer', 'inducer_conc', 'temperature', 'perturbation_idx']
model_params = ['alpha_rep', 'm_peri', 'phi_peri', 'rho_peri']
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
                                   'inducer_conc', 'inducer', 'temperature', 'perturbation_idx']):
        percs.loc[percs[f'{p}_dim_0'] == g[-1]-1, ['strain', 'carbon_source',
                                                   'overexpression', 'inducer_conc', 'inducer', 'temperature']] = g[:-1]
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
                                   'inducer', 'inducer_conc', 'temperature', 'perturbation_idx']):
        percs.loc[percs[f'{p}_dim_0'] == g[-1]-1, ['strain', 'carbon_source', 'overexpression',
                                                   'inducer', 'inducer_conc', 'temperature']] = g[:-1]
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
for g, d in growth_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer', 'inducer_conc', 'temperature', 'perturbation_idx']):
    percs.loc[percs[f'growth_mu_dim_0'] == g[-1]-1, ['strain', 'carbon_source', 'overexpression',
                                                     'inducer', 'inducer_conc', 'temperature']] = g[:-1]
percs.drop(columns=f'growth_mu_dim_0', inplace=True)
percs.to_csv('../../data/mcmc/growth_parameter_percentiles.csv', index=False)
