# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
import size.viz
import seaborn as sns
import scipy.stats
cor, pal = size.viz.matplotlib_style()

# Compile the statistical model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/one-shot_complete_inference.stan')

# %%
# Load the various datasets
cal_data = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
brad_data = pd.read_csv(
    '../../data/protein_quantification/bradford_periplasmic_protein.csv')
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
biomass_data = pd.read_csv(
    '../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')

# Keep only the bradford data with more than two replicates
brad_data = pd.concat([d for _, d in brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']) if len(d) > 2], sort=False)
size_data = size_data[size_data['temperature_C'] == 37]
size_data = pd.concat([d for _, d in size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']) if len(d) > 2], sort=False)

# %%
# Filter, label, and transform bradford data
# brad_data = brad_data[brad_data['strain'].isin(
#     ['wildtype', 'malE-rbsB-fliC-KO'])]
brad_data['cond_idx'] = brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']).ngroup() + 1
brad_data['conv_factor'] = brad_data['dilution_factor'] * \
    brad_data['extraction_volume_mL'] / \
    (brad_data['od_600nm'] * brad_data['culture_volume_mL'])
brad_data['od_per_biomass'] = brad_data['od_595nm'] * brad_data['conv_factor']

# Add indexing to the size data
size_data['size_cond_idx'] = size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']).ngroup() + 1


# Define teh data dictionary
data_dict = {'N_cal': len(cal_data),
             'concentration': cal_data['protein_conc_ug_ml'].values.astype(float),
             'cal_od': cal_data['od_595nm'].values.astype(float),

             'N_brad': len(brad_data),
             'J_brad_cond': brad_data['cond_idx'].max(),
             'brad_cond_idx': brad_data['cond_idx'].values.astype(int),
             'brad_od': brad_data['od_595nm'].values.astype(float),
             'conv_factor': brad_data['conv_factor'].values.astype(float),

             'N_size': len(size_data),
             'J_size_cond': size_data['size_cond_idx'].max(),
             'size_cond_idx': size_data['size_cond_idx'].values.astype(int),
             'width': size_data['width_median'].values.astype(float),
             'length': size_data['length'].values.astype(float),
             'volume': size_data['volume'].values.astype(float),
             'peri_volume': size_data['periplasm_volume'].values.astype(float),
             'surface_area': size_data['surface_area'].values.astype(float),
             'surface_area_volume': size_data['surface_to_volume'].values.astype(float),

             'N_biomass': len(biomass_data),
             'biomass': biomass_data['dry_mass_fg']
             }


# %%
# Sample the posterior
_samples = model.sample(data_dict)
samples = az.from_cmdstanpy(_samples)

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
ax.set_xlabel('OD$_{595nm}$ / OD$_{600nm} \cdot$ mL')

# %%
# Biomass ppc
biomass_ppc = samples.posterior.biomass_rep.to_dataframe()
biomass_ppc['idx'] = 1
biomass_percs = size.viz.compute_percentiles(biomass_ppc, 'biomass_rep', 'idx')

fig, ax = plt.subplots(1, 1, figsize=(3, 1))

for i, (g, d) in enumerate(biomass_percs.groupby(['interval'], sort=False)):
    ax.hlines(1, d['lower'], d['upper'], lw=15, zorder=i+1, color=ppc_cmap[g])


ax.plot(biomass_data['dry_mass_fg'], np.ones(
    len(biomass_data)), 'o', color=cor['primary_green'], zorder=i+1)
ax.set_yticks([])
ax.set_xlabel('dry mass [fg / OD$_{600nm}$ mL]')

# %%
fig, ax = plt.subplots(2, 3, figsize=(6, 12), sharey=True)
props = ['width', 'length', 'volume', 'peri_volume',
         'surface_area', 'surface_area_vol']
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
labels = [g[1:] for g, _ in size_data.groupby(
    ['size_cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])]
_ = ax[0, 0].set_yticks(size_data['size_cond_idx'].unique())
_ = ax[1, 0].set_yticks(size_data['size_cond_idx'].unique())
_ = ax[0, 0].set_yticklabels(labels)
_ = ax[1, 0].set_yticklabels(labels)


# %%
# Link the axes
brad_idx = {g[1:]: g[0] for g, _ in brad_data.groupby(
    ['cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL'])}
size_loc = {g[0]: brad_idx[g[1:]] for g, _ in size_data.groupby(
    ['size_cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc']) if g[1:] in brad_idx.keys()}
width_ppc = samples.posterior.width_mu.to_dataframe().reset_index()

frac_ppc = samples.posterior.periplasmic_biomass_fraction.to_dataframe().reset_index()

for g, _ in brad_data.groupby(['cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']):
    frac_ppc.loc[frac_ppc['periplasmic_biomass_fraction_dim_0'] == g[0]-1, [
        'strain', 'carbon_source', 'overexpression', 'inducer_conc']] = g[1:]
    frac_ppc.loc[frac_ppc['periplasmic_biomass_fraction_dim_0']
                 == g[0]-1, 'link_idx'] = brad_idx[g[1:]]
frac_perc = size.viz.compute_percentiles(frac_ppc, 'periplasmic_biomass_fraction', [
                                         'link_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])

for g, _ in size_data.groupby(['size_cond_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc']):
    if g[0] in size_loc:
        width_ppc.loc[width_ppc['width_mu_dim_0'] == g[0]-1, [
            'strain', 'carbon_source', 'overexpression', 'inducer_conc']] = g[1:]
        width_ppc.loc[width_ppc['width_mu_dim_0']
                      == g[0]-1, 'link_idx'] = size_loc[g[0]]
width_perc = size.viz.compute_percentiles(width_ppc, 'width_mu', [
                                          'link_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])
width_ppc.dropna(inplace=True)

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
