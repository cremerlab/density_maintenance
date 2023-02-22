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

# Load the necessary datasets
bradford_cal_data = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
bradford_prot_data = pd.read_csv(
    '../../data/protein_quantification/bradford_periplasmic_protein.csv')
growth_data = pd.read_csv(
    '../../data/growth_curves/growth_measurements_processed.csv')
flow_data = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
lit_data = pd.read_csv(
    '../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')

flow_data = flow_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc', 'date', 'run_no'])[
    'cells_per_biomass'].mean().reset_index()

flow_data = flow_data[flow_data['strain'] == 'wildtype']

# Restrict to temperature = 37 and drop the minKO and d2 samples
size_data = size_data[size_data['temperature_C'] == 37]
size_data = size_data[size_data['strain'].isin(['wildtype', 'malE-rbsB-fliC-KO',
                                                'lpp14', 'lpp21'])].copy()
# Drop size measurements with only one replicate
# counted = size_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc']).count().reset_index()
include = [d for _, d in size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']) if len(d) > 1]
size_data = pd.concat(include, sort=False)

bradford_prot_data = bradford_prot_data[bradford_prot_data['strain'].isin(
    ['wildtype', 'malE-rbsB-fliC-KO', 'lpp14', 'lpp21'])].copy()
bradford_prot_data['conv_factor'] = bradford_prot_data['dilution_factor'] * \
    bradford_prot_data['extraction_volume_mL'] / \
    (bradford_prot_data['od_600nm'] * bradford_prot_data['culture_volume_mL'])
# Add appropriate indexing
size_data['idx'] = size_data.groupby(['strain', 'carbon_source',
                                      'overexpression', 'inducer_conc']).ngroup() + 1

# Set up a dictionary mapping idx to details
idx_mapper = {}
for g, _ in size_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc', 'idx']):
    idx_mapper[tuple(g[:-1])] = g[-1]

# Add mapping to protein quantification data
bradford_prot_data['link_idx'] = 0
for g, d in bradford_prot_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']):
    if tuple(g) in idx_mapper.keys():
        bradford_prot_data.loc[(bradford_prot_data['strain'] == g[0]) &
                               (bradford_prot_data['carbon_source'] == g[1]) &
                               (bradford_prot_data['overexpression'] == g[2]) &
                               (bradford_prot_data['inducer_conc_ng_mL'] == g[3]),
                               'link_idx'] = idx_mapper[tuple(g)]

# Drop protein data for which there are no size measurements
bradford_prot_data = bradford_prot_data[bradford_prot_data['link_idx'] != 0].copy(
)

# Add mapping to flow cytometry data
flow_data['link_idx'] = 0
for g, d in flow_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc']):
    if tuple(g) in idx_mapper.keys():
        flow_data.loc[(flow_data['strain'] == g[0]) &
                      (flow_data['carbon_source'] == g[1]) &
                      (flow_data['overexpression'] == g[2]) &
                      (flow_data['inducer_conc'] == g[3]),
                      'link_idx'] = idx_mapper[tuple(g)]

# Add mapping to growth rate data
growth_data['link_idx'] = 0
for g, d in growth_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_ml']):
    if tuple(g) in idx_mapper.keys():
        growth_data.loc[(growth_data['strain'] == g[0]) &
                        (growth_data['carbon_source'] == g[1]) &
                        (growth_data['overexpression'] == g[2]) &
                        (growth_data['inducer_conc_ng_ml'] == g[3]),
                        'link_idx'] = idx_mapper[tuple(g)]

# Drop growth curves for which we don't have size data
growth_data = growth_data[growth_data['link_idx'] != 0].copy()


# %%
# Add intrasample idx
bradford_prot_data['brad_cond_idx'] = bradford_prot_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']).ngroup() + 1
flow_data['flow_cond_idx'] = flow_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']).ngroup() + 1
growth_data['growth_cond_idx'] = growth_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_ml']).ngroup() + 1
growth_data['growth_curve_idx'] = growth_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_ml',
     'run_idx']).ngroup() + 1


# %%
# Generate the data dictionary
data_dict = {
    'J_size_cond': size_data['idx'].max(),
    'J_growth_cond': growth_data['growth_cond_idx'].max(),
    'J_growth_curves': growth_data['growth_curve_idx'].max(),
    'J_brad_cond': bradford_prot_data['brad_cond_idx'].max(),
    'J_flow_cond': flow_data['flow_cond_idx'].max(),

    'N_brad': len(bradford_prot_data),
    'N_flow': len(flow_data),
    'N_growth': len(growth_data),
    'N_size': len(size_data),

    'brad_link_idx': bradford_prot_data['link_idx'].values.astype(int),
    'flow_link_idx': flow_data['link_idx'].values.astype(int),
    'growth_link_idx': growth_data.groupby(['growth_curve_idx'])['link_idx'].min().astype(int),
    'size_idx': size_data['idx'].values.astype(int),

    'brad_idx': bradford_prot_data['brad_cond_idx'].values.astype(int),
    'flow_idx': flow_data['flow_cond_idx'].values.astype(int),
    'growth_cond_idx': growth_data.groupby(['growth_curve_idx'])['growth_cond_idx'].min().astype(int),
    'growth_curve_idx': growth_data['growth_curve_idx'].values.astype(int),

    'N_cal_meas': len(bradford_cal_data),
    'concentration': bradford_cal_data['protein_conc_ug_ml'].values.astype(float),
    'cal_od': bradford_cal_data['od_595nm'].values.astype(float),

    'prot_od': bradford_prot_data['od_595nm'].values.astype(float),
    'od_conv_factor': bradford_prot_data['conv_factor'].values.astype(float),

    'N_lit_meas': len(lit_data),
    'drymass': lit_data['dry_mass_fg'].values.astype(float),

    'growth_time': growth_data['elapsed_time_hr'].values.astype(float),
    'growth_od': growth_data['od_600nm'].values.astype(float),

    'cells_per_biomass': flow_data['cells_per_biomass'].values.astype(float),

    'mean_width': size_data['width_median'].values.astype(float),
    'mean_length': size_data['length'].values.astype(float),
    'mean_vol': size_data['volume'].values.astype(float),
    'mean_peri_vol': size_data['periplasm_volume'].values.astype(float),
    'mean_peri_vol_frac': size_data['periplasm_volume_fraction'].values.astype(float),
    'mean_sa': size_data['surface_area'].values.astype(float),
    'mean_sav': size_data['surface_to_volume'].values.astype(float),
    'mean_width_mean': size_data.groupby(['idx'])['width_median'].mean().astype(float),
    'mean_width_std': size_data.groupby(['idx'])['width_median'].std().astype(float),
    'mean_length_mean': size_data.groupby(['idx'])['length'].mean().astype(float),
    'mean_length_std': size_data.groupby(['idx'])['length'].std().astype(float),
    'mean_vol_mean': size_data.groupby(['idx'])['volume'].mean().astype(float),
    'mean_vol_std': size_data.groupby(['idx'])['volume'].std().astype(float),
    'mean_peri_vol_mean': size_data.groupby(['idx'])['periplasm_volume'].mean().astype(float),
    'mean_peri_vol_std': size_data.groupby(['idx'])['periplasm_volume'].std().astype(float),
    'mean_peri_vol_frac_mean': size_data.groupby(['idx'])['periplasm_volume_fraction'].mean().astype(float),
    'mean_peri_vol_frac_std': size_data.groupby(['idx'])['periplasm_volume_fraction'].std().astype(float),
    'mean_sa_mean': size_data.groupby(['idx'])['surface_area'].mean().astype(float),
    'mean_sa_std': size_data.groupby(['idx'])['surface_area'].std().astype(float),
    'mean_sav_mean': size_data.groupby(['idx'])['surface_to_volume'].mean().astype(float),
    'mean_sav_std': size_data.groupby(['idx'])['surface_to_volume'].std().astype(float),
}

# %%
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/one-shot_semicomplete_inference.stan')

# %%
_samples = model.sample(data_dict)
samples = az.from_cmdstanpy(_samples)


# %%
# Calibration curve ppc
def compute_percentiles(df,
                        quantity,
                        groupby,
                        lower_bounds=[5, 10, 15, 20, 25, 30, 35, 40, 45, ],
                        upper_bounds=[95, 90, 85, 80, 75, 70, 65, 60, 55],
                        interval_labels=['90%', '80%', '70%', '60%',
                                         '50%', '40%', '30%', '20%', '10%']):

    # Allow flexibility in what quantities are being supplied
    if type(quantity) != str:
        if type(quantity) != list:
            raise TypeError("`quantity` must be a `str` or list of `str.`")
    else:
        quantity = [quantity]

    # Instantiate the dataframe and loop through every group
    perc_df = pd.DataFrame([])
    if type(groupby) != list:
        groupby = [groupby]

    for g, d in df.groupby(groupby):
        if type(g) != tuple:
            g = (g,)
        # Compute the percentiles for different quantities
        for q in quantity:
            lower = np.percentile(d[f'{q}'].values, lower_bounds)
            upper = np.percentile(d[f'{q}'].values, upper_bounds)
            _df = pd.DataFrame(np.array([lower, upper]).T, columns=[
                               'lower', 'upper'])
            _df['quantity'] = q
            _df['interval'] = interval_labels

            # Add the grouping informaton
            for i in range(len(groupby)):
                _df[groupby[i]] = g[i]
            perc_df = pd.concat([perc_df, _df], sort=False)
    return perc_df


ppc_percs = ['90%', '80%', '70%', '60%', '50%', '40%', '30%', '20%', '10%']
ppc_cmap = sns.color_palette('Greys_r', n_colors=len(ppc_percs) + 2).as_hex()
ppc_cmap_dict = {i: c for i, c in zip(
    [f'{k}%' for k in np.arange(10, 100, 10)], ppc_cmap)}

# # %%
# # Unpack the calbiration curve ppc
concs = bradford_cal_data['protein_conc_ug_ml'].values
cal_ppc_df = samples.posterior.od595_calib_rep.to_dataframe().reset_index()
for i, c in enumerate(concs):
    cal_ppc_df.loc[cal_ppc_df['od595_calib_rep_dim_0']
                   == i, 'protein_conc_ug_ml'] = c


# PPC for calibration curve
cal_ppc_percs = compute_percentiles(
    cal_ppc_df, 'od595_calib_rep', 'protein_conc_ug_ml')
cal_ppc_percs

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
for i, (g, d) in enumerate(cal_ppc_percs.groupby('interval', sort=False)):
    ax.fill_between(d['protein_conc_ug_ml'], d['lower'], d['upper'],
                    color=ppc_cmap_dict[g], zorder=i+1, label=g)

ax.plot(bradford_cal_data['protein_conc_ug_ml'], bradford_cal_data['od_595nm'], 'o', color=cor['primary_red'],
        label='measurement', ms=3, zorder=1000)
ax.legend()
ax.set_xlabel('protein concentration [µg / mL]')
ax.set_ylabel('OD$_{595nm}$ [a.u.]')
ax.set_title('Bradford assay calibration curve')
plt.tight_layout()
# plt.savefig('../../figures/mcmc/one-shot_calibration_curve_ppc.pdf')
# %%

# # PPC for bradford assay
# prot_ppc_df = samples.posterior.od595_per_biomass_rep.to_dataframe().reset_index()
# # Map the index to the carbon source
# _prot_data = bradford_prot_data[['carbon_source', 'cond_idx']]
# _prot_data['n'] = np.arange(len(_prot_data))
# for k, v in zip(_prot_data['carbon_source'].values, _prot_data['n'].values):
#     prot_ppc_df.loc[prot_ppc_df['od595_per_biomass_rep_dim_0']
#                     == v, 'carbon_source'] = k

# prot_ppc_df

# Compute the percentiles
prot_ppc_percs = compute_percentiles(
    prot_ppc_df, 'od595_per_biomass_rep', 'carbon_source')

# Plot the PPC under the measurements
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
for i, (g, d) in enumerate(prot_ppc_percs.groupby(['carbon_source', 'interval'], sort=False)):
    ax.hlines(carb_loc[g[0]], d['lower'], d['upper'],
              color=ppc_cmap_dict[g[1]], lw=15, zorder=i)

for g, d in bradford_prot_data.groupby(['carbon_source']):
    ax.plot(d['od_595nm'] * d['conv_factor'], carb_loc[g] + np.random.normal(0, 0.01,
            len(d)), 'o', ms=3, color=cor['primary_red'], zorder=1000)

_ = ax.set_yticks(list(carb_loc.values())[:-1])
_ = ax.set_yticklabels(list(carb_loc.keys())[:-1])
_ = ax.set_xlabel('OD$_{595nm}$ per biomass')
ax.set_title('Periplasmic protein quantification measurements')
plt.tight_layout()
plt.savefig('../../figures/mcmc/one-shot_periprot_quantification_ppc.pdf')


# # %%
# # Lit data ppc
# lit_data_ppc = samples.posterior.drymass_rep.to_dataframe().reset_index()
# lit_data_ppc['idx'] = 1
# lit_data_percs = compute_percentiles(lit_data_ppc, 'drymass_rep', 'idx')

# fig, ax = plt.subplots(1, 1, figsize=(4, 1))
# for i, (g, d) in enumerate(lit_data_percs.groupby(['interval'], sort=False)):
#     ax.hlines(0, d['lower'], d['upper'], zorder=i +
#               1, color=ppc_cmap_dict[g], lw=15)

# ax.plot(lit_data['dry_mass_fg'], np.ones(len(lit_data))*0, 'o', ms=4,
#         color=cor['primary_green'], zorder=1000)

# ax.set_yticks([])
# ax.set_xlabel('total drymass [µg/OD$_{600nm} \cdot$mL')
# ax.set_title('Basan et al. 2015 drymass measurements')
# plt.tight_layout()
# plt.savefig('../../figures/mcmc/one-shot_drymass_quantification_ppc.pdf')

# %%
cpb_ppc_df = samples.posterior[[
    'cells_per_biomass_rep']].to_dataframe().reset_index()

# Map conditions to carbon_source
flow_data['n'] = np.arange(len(flow_data))
for dim, carb, strain in zip(flow_data['n'].values, flow_data['carbon_source'].values, flow_data['strain'].values):
    cpb_ppc_df.loc[cpb_ppc_df['cells_per_biomass_rep_dim_0']
                   == dim, 'carbon_source'] = carb
    cpb_ppc_df.loc[cpb_ppc_df['cells_per_biomass_rep_dim_0']
                   == dim, 'strain'] = strain

# # Compute the percentiles
cpb_perc_df = compute_percentiles(
    cpb_ppc_df, 'cells_per_biomass_rep', ['carbon_source', 'strain'])
fig, ax = plt.subplots(1, 1, figsize=(4, 2))

carb_loc = {k: i for i, k in enumerate(
    ['LB', 'glucose', 'glucoseCAA', 'sorbitol', 'glycerol', 'acetate'])}
for i, (g, d) in enumerate(cpb_perc_df[cpb_perc_df['strain'] == 'wildtype'].groupby(['carbon_source', 'strain', 'interval'], sort=False)):
    # ax.plot(d['lower'], carb_loc[g[0]], 'o', color=ppc_cmap_dict[g[2]])
    # ax.plot(d['upper'], carb_loc[g[0]], 'o', color=ppc_cmap_dict[g[2]])
    ax.hlines(carb_loc[g[0]], d['lower'], d['upper'],
              lw=10, color=ppc_cmap_dict[g[2]])
# ax.set_xscale('log')

for g, d in flow_data[flow_data['strain'] == 'wildtype'].groupby(['carbon_source']):
    ax.plot(d['cells_per_biomass'], carb_loc[g] + np.random.normal(0, 0.1, len(d)), 'o',
            color=cor['primary_red'], ms=4, zorder=1000)
# %%
# # ax.set_xlim([1E8, 5E9])
# ax.set_yticks([1, 2, 3, 4, 5, 6])
# ax.set_yticklabels(['acetate', 'sorbitol', 'glycerol',
#                    'glucose', 'glucoseCAA', 'LB'])
# ax.set_xlabel('cells per biomass')
# ax.set_title('Flow cytometry event counts')
# plt.tight_layout()
# plt.savefig('../../figures/mcmc/one-shot_flow_cytometry_ppc.pdf')


# # %%
# # Plot the fit to the flow data
# cpb_fit_df = samples.posterior[[
#     'vol_cells_per_biomass']].to_dataframe().reset_index()
# for carb, dim in carb_idx_dict.items():
#     cpb_fit_df.loc[cpb_fit_df['vol_cells_per_biomass_dim_0']
#                    == dim-1, 'carbon_source'] = carb

# _df = samples.posterior[[
#     'vol_mu']].to_dataframe().reset_index()
# cpb_fit_df['vol_mu'] = _df['vol_mu'].values

# # Compute the percentiles
# cpb_fit_percs = compute_percentiles(
#     cpb_fit_df, ['vol_cells_per_biomass', 'vol_mu'], 'carbon_source')
# cpb_fit_percs

# fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# for g, d in cpb_fit_percs.groupby(['carbon_source'], sort=False):
#     lam = d[d['quantity'] == 'vol_mu']
#     cpb = d[d['quantity'] == 'vol_cells_per_biomass']
#     perc_10_cpb = cpb[cpb['interval'] == '10%']
#     perc_10_lam = lam[lam['interval'] == '10%']
#     x = np.mean([perc_10_lam['lower'].values[0],
#                 perc_10_lam['upper'].values[0]])
#     y = np.mean([perc_10_cpb['lower'].values[0],
#                 perc_10_cpb['upper'].values[0]])
#     for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
#         _lam = _d[_d['quantity'] == 'vol_mu']
#         _cpb = _d[_d['quantity'] == 'vol_cells_per_biomass']
#         ax.vlines(x, _cpb['lower'], _cpb['upper'], lw=3,
#                   zorder=i+1, color=ppc_cmap_dict[_g])
#         ax.hlines(y, _lam['lower'], _lam['upper'], lw=3,
#                   zorder=i+1, color=ppc_cmap_dict[_g])

#     if g in flow_data['carbon_source'].values:
#         _data = flow_data[flow_data['carbon_source'] == g]
#         ax.plot(x * np.ones(len(_data)),
#                 _data['cells_per_biomass'], 'o', ms=5, color=cor['primary_red'], zorder=1000)

# # compute and plot the main fit
# fit_params = samples.posterior[[
#     'k_cells_per_biomass_tilde', 'beta_0_tilde']].to_dataframe().reset_index()
# vol_range = np.linspace(0.1, 4, 200)
# fit = 1E9 * (fit_params['beta_0_tilde'].mean() -
#              fit_params['k_cells_per_biomass_tilde'].mean() * vol_range)
# ax.plot(vol_range, fit, '-', lw=1, color=cor['primary_blue'])
# ax.set_xlabel('average cell volume [µm$^{3}$]')
# ax.set_ylabel('cells per biomass')
# ax.set_title('Linear volume-dependence on cell density')

# plt.tight_layout()
# plt.savefig('../../figures/mcmc/one-shot_cell_density_vol_dependence.pdf')

# # %%
# # Plot the growth ppcs
# lam_ppc_df = samples.posterior.growth_od_rep.to_dataframe().reset_index()

# # Map time and carbon sources to everything
# carbs = growth_data['carbon_source'].values
# breps = growth_data['brep_idx'].values
# elapsed_time = growth_data['elapsed_time_hr'].values
# dims = np.arange(len(growth_data))
# for dim, c, b, t in zip(dims, carbs, breps, elapsed_time):
#     lam_ppc_df.loc[lam_ppc_df['growth_od_rep_dim_0']
#                    == dim, 'carbon_source'] = c
#     lam_ppc_df.loc[lam_ppc_df['growth_od_rep_dim_0'] == dim, 'brep'] = int(b)
#     lam_ppc_df.loc[lam_ppc_df['growth_od_rep_dim_0'] == dim, 'time'] = t
# lam_ppc_percs = compute_percentiles(lam_ppc_df, 'growth_od_rep', [
#                                     'carbon_source', 'brep', 'time'])
# # %%
# for g, d in growth_data.groupby(['carbon_source']):
#     n_reps = len(d[d['elapsed_time_hr'] == 0])
#     cols = 3
#     rows = int(np.ceil(n_reps / cols))
#     ex = cols * rows - n_reps
#     fig, ax = plt.subplots(rows, cols, figsize=(8, 2 * rows))

#     # Turn off unpopulated axes
#     if ex != 0:
#         for a in ax.ravel()[-ex:]:
#             a.axis(False)
#         ax = ax.ravel()[:-ex]
#     else:
#         ax = ax.ravel()
#     for a in ax:
#         a.set_xlabel('elapsed time [hr]')
#         a.set_ylabel('log OD$_{600nm}$')

#     # get the ppc
#     _ppc = lam_ppc_percs[(lam_ppc_percs['carbon_source'] == g)]
#     ind = 0
#     for i, (_g, _d) in enumerate(_ppc.groupby(['brep', 'interval'], sort=False)):
#         if i == 0:
#             _brep = _g[0]
#         if _g[0] != _brep:
#             ind += 1
#             _brep = _g[0]

#         ax[ind].fill_between(_d['time'], np.log(_d['lower']), np.log(_d['upper']),
#                              color=ppc_cmap_dict[_g[1]], zorder=i+1)

#     for i, (_g, _d) in enumerate(d.groupby(['brep_idx'])):
#         ax[i].plot(_d['elapsed_time_hr'], np.log(_d['od_600nm']),
#                    '-o', color=cor['primary_red'], ms=4, zorder=1000, lw=1)
#     plt.suptitle(f'wildtype {g} growth curves')
#     plt.tight_layout()
#     plt.savefig(f'../../figures/mcmc/one-shot_{g}_growth_curves_ppc.pdf')
#     plt.close()

# # %%
# Cell Size ppc -- need to do this piecewise
size_perc_df = pd.DataFrame([])
for val in ['width_rep', 'length_rep', 'vol_rep', 'peri_vol_rep',
            'peri_vol_frac_rep', 'sa_rep', 'sav_rep']:
    _ppc_df = samples.posterior[f'{val}'].to_dataframe().reset_index()
    # Map the carbon sources
    for dim, carb in zip(np.arange(len(size_data)), size_data['carbon_source'].values):
        _ppc_df.loc[_ppc_df[f'{val}_dim_0'] == dim, 'carbon_source'] = carb

    # Compute percentiles
    perc_df = compute_percentiles(_ppc_df, f'{val}', ['carbon_source'])
    size_perc_df = pd.concat([size_perc_df, perc_df], sort=False)

# %%
fig, ax = plt.subplots(2, 4, figsize=(8, 4), sharey=True)
labs = ['width [µm]', 'length [µm]', 'volume [µm$^3$]', 'periplasm volume\n[µm$^3$]',
        'periplasm volume fraction', 'surface area [µm$^2$]', 'surface-to-volume [µm$^{-1}$] ']
ax = ax.ravel()
for i, ell in enumerate(labs):
    ax[i].set_xlabel(ell)

for a in [ax[0], ax[4]]:
    a.set_yticks(list(carb_loc.values()))
    a.set_yticklabels(list(carb_loc.keys()))


# Plot the percentiles
for i, (g, d) in enumerate(size_perc_df.groupby(['quantity'], sort=False)):
    for j, (_g, _d) in enumerate(d.groupby(['carbon_source', 'interval'], sort=False)):
        ax[i].hlines(carb_loc[_g[0]], _d['lower'], _d['upper'],
                     lw=15, color=ppc_cmap_dict[_g[1]])

# Plot the data
for g, d in size_data.groupby(['carbon_source']):
    for i, v in enumerate(['width_median', 'length', 'volume', 'periplasm_volume',
                           'periplasm_volume_fraction', 'surface_area', 'surface_to_volume']):
        ax[i].plot(d[f'{v}'], carb_loc[g] * np.ones(len(d)) +
                   np.random.normal(0, 0.05, len(d)), 'o', ms=4, color=cor['primary_red'])
ax[-1].axis(False)
plt.suptitle('average cell dimension measurements',
             fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('../../figures/mcmc/one-shot_cell_dimension_ppc.pdf')
# # %%
# sav_rho_percs = pd.DataFrame([])
# for i, k in enumerate(['peri_drymass_frac', 'sav_mu']):
#     post_df = samples.posterior[f'{k}'].to_dataframe().reset_index()
#     for key, val in carb_idx_dict.items():
#         post_df.loc[post_df[f'{k}_dim_0'] == val - 1, 'carbon_source'] = key
#     perc = compute_percentiles(post_df, k, 'carbon_source')
#     perc = perc[perc['carbon_source'] != 'LB']
#     sav_rho_percs = pd.concat([sav_rho_percs, perc], sort=False)


# fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# for g, d in sav_rho_percs.groupby(['carbon_source']):
#     if g == 'LB':
#         continue
#     sav_10 = d[(d['quantity'] == 'sav_mu') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     frac_10 = d[(d['quantity'] == 'peri_drymass_frac') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
#         sav = _d[_d['quantity'] == 'sav_mu']
#         frac = _d[_d['quantity'] == 'peri_drymass_frac']
#         ax.vlines(frac_10 * 100, sav['lower'], sav['upper'],
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+1)
#         ax.hlines(sav_10, frac['lower'] * 100, frac['upper'] * 100,
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+1)


# # Plot the theory line
# k = (0.03/0.15) * (1 - 0.15)/(1 - 0.03)
# delta = 0.025
# phi_range = np.linspace(0, 0.05, 200)
# phi_term = k * (1 - phi_range) / phi_range
# theory = (delta * (phi_term))**-1
# ax.plot(phi_range * 100, theory, 'k--')
# # ax.set_ylim([4, 8])
# ax.set_xlabel('periplasmic protein drymass fraction [%]')
# ax.set_ylabel('surface-to-volume [µm$^{-1}$]')

# # %%
mass_fracs = pd.read_csv('../../data/compiled_mass_fractions.csv')
genes = pd.read_csv('./genes_classification_all.csv')
_genes = genes[genes['location'].isin(['IM', 'OM', 'PE', 'LPO'])]
mass_fracs = mass_fracs[mass_fracs['gene_name'].isin(_genes['gene'].unique())]

mass_fracs
# # %%
# gr = samples.posterior.growth_mu.to_dataframe().reset_index()
# for carb, idx in carb_idx_dict.items():
#     gr.loc[gr['growth_mu_dim_0'] == idx-1, 'carbon_source'] = carb
# gr_agg = gr.groupby(['carbon_source'])['growth_mu'].median().reset_index()

# for carb, lam in zip(gr_agg['carbon_source'].values, gr_agg['growth_mu'].values):
#     size_data.loc[size_data['carbon_source'] == carb, 'growth_mu'] = lam
# # %%
# peri_mass_fracs = mass_fracs[mass_fracs['periplasm'] == True].groupby(
#     ['dataset_name', 'condition', 'growth_rate_hr'])['mass_frac'].sum().reset_index()

# dry_frac = 0.3
# prot_frac = 0.55
# density = 1.1

# # COmpute simple fits
w_popt = scipy.stats.linregress(
    size_data['growth_mu'], size_data['width_median'])
ell_popt = scipy.stats.linregress(
    size_data['growth_mu'], np.log(size_data['length']))
peri_vol_popt = scipy.stats.linregress(
    size_data['growth_mu'], np.log(size_data['periplasm_volume']))
vol_popt = scipy.stats.linregress(
    size_data['growth_mu'], np.log(size_data['volume']))
sav_popt = scipy.stats.linregress(
    size_data['growth_mu'], np.log(size_data['surface_to_volume']))

# Compute the periplasmic protein density
peri_mass_fracs['width'] = w_popt[0] * \
    peri_mass_fracs['growth_rate_hr'] + w_popt[1]
peri_mass_fracs['length'] = np.exp(
    ell_popt[0] * peri_mass_fracs['growth_rate_hr'] + ell_popt[1])
peri_mass_fracs['peri_vol'] = np.exp(
    peri_vol_popt[0] * peri_mass_fracs['growth_rate_hr'] + peri_vol_popt[1])
peri_mass_fracs['volume'] = np.exp(
    vol_popt[0] * peri_mass_fracs['growth_rate_hr'] + vol_popt[1])
peri_mass_fracs['sav'] = np.exp(
    sav_popt[0] * peri_mass_fracs['growth_rate_hr'] + sav_popt[1])
# peri_mass_fracs['peri_volume'] = size.analytical.surface_area(peri_mass_fracs['length'], peri_mass_fracs['width']) * 0.025
peri_mass_fracs['tot_protein'] = density * \
    dry_frac * prot_frac * peri_mass_fracs['volume']
peri_mass_fracs['peri_protein'] = peri_mass_fracs['mass_frac'] * \
    peri_mass_fracs['tot_protein']
peri_mass_fracs['rho_peri'] = (
    peri_mass_fracs['peri_protein'] * 1E3) / peri_mass_fracs['peri_vol']
peri_mass_fracs['biomass_frac'] = peri_mass_fracs['peri_protein'] / \
    (

# # %%

# # %%
# fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# ax.plot(peri_mass_fracs['biomass_frac']*100,
#         peri_mass_fracs['sav'], 'X', alpha=0.2)

# for g, d in sav_rho_percs.groupby(['carbon_source']):
#     if g == 'LB':
#         continue
#     sav_10 = d[(d['quantity'] == 'sav_mu') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     frac_10 = d[(d['quantity'] == 'peri_drymass_frac') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
#         sav = _d[_d['quantity'] == 'sav_mu']
#         frac = _d[_d['quantity'] == 'peri_drymass_frac']
#         ax.vlines(frac_10 * 100, sav['lower'], sav['upper'],
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+8)
#         ax.hlines(sav_10, frac['lower'] * 100, frac['upper'] * 100,
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+8)


# # Plot the theory line
# k = (0.03/0.15) * (1 - 0.15)/(1 - 0.03)
# delta = 0.025
# phi_range = np.linspace(0, 0.08, 200)
# phi_term = k * (1 - phi_range) / phi_range
# theory = (delta * (phi_term))**-1
# ax.plot(phi_range * 100, theory, 'k--')
# # ax.set_ylim([4, 8])
# ax.set_xlim([0, 8])
# ax.set_xlabel('periplasmic protein drymass fraction [%]')
# ax.set_ylabel('surface-to-volume [µm$^{-1}$]')
# plt.savefig('../../figures/mcmc/one-shot_theory_comparison.pdf')

# # %%
# rho_percs = pd.DataFrame([])
# for i, k in enumerate(['peri_density', 'growth_mu']):
#     post_df = samples.posterior[f'{k}'].to_dataframe().reset_index()
#     for key, val in carb_idx_dict.items():
#         post_df.loc[post_df[f'{k}_dim_0'] == val - 1, 'carbon_source'] = key
#     perc = compute_percentiles(post_df, k, 'carbon_source')
#     rho_percs = pd.concat([rho_percs, perc], sort=False)

#  # %%
# fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# ax.plot(peri_mass_fracs['growth_rate_hr'],
#         peri_mass_fracs['rho_peri'], 'X', alpha=0.5)
# for g, d in rho_percs.groupby(['carbon_source']):
#     if g == 'LB':
#         continue
#     lam_10 = d[(d['quantity'] == 'growth_mu') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     rho_10 = d[(d['quantity'] == 'peri_density') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
#         lam = _d[_d['quantity'] == 'growth_mu']
#         rho = _d[_d['quantity'] == 'peri_density']
#         ax.hlines(rho_10, lam['lower'], lam['upper'],
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+1)
#         ax.vlines(lam_10, rho['lower'], rho['upper'],
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+1)

# ax.set_xlabel('growth rate [hr$^{-1}$]')
# ax.set_ylabel('periplasmic protein density [fg/fL]')
# ax.set_ylim([0, 175])
# ax.set_title('Comparison with mass spectrometry data')
# # plt.savefig('../../figures/mcmc/one-shot_peri_density_comparison.pdf')

# # %%
# ppc_cmap = sns.color_palette('Reds_r', n_colors=len(
#     cal_ppc_percs['interval'].unique()) + 2).as_hex()
# ppc_cmap_dict = {i: c for i, c in zip(
#     [f'{k}%' for k in np.arange(10, 100, 10)], ppc_cmap)}

# width_rho_percs = pd.DataFrame([])
# for i, k in enumerate(['peri_drymass_frac', 'width_mu']):
#     post_df = samples.posterior[f'{k}'].to_dataframe().reset_index()
#     for key, val in carb_idx_dict.items():
#         post_df.loc[post_df[f'{k}_dim_0'] == val - 1, 'carbon_source'] = key
#     perc = compute_percentiles(post_df, k, 'carbon_source')
#     perc = perc[perc['carbon_source'] != 'LB']
#     width_rho_percs = pd.concat([width_rho_percs, perc], sort=False)


# fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# ax.plot(peri_mass_fracs['biomass_frac']*100,
#         1/peri_mass_fracs['width'], 'X', alpha=0.2, color=cor['primary_black'], zorder=1, ms=4, label='__nolegend__')
# for g, d in width_rho_percs.groupby(['carbon_source']):
#     if g == 'LB':
#         continue
#     width_10 = d[(d['quantity'] == 'width_mu') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     frac_10 = d[(d['quantity'] == 'peri_drymass_frac') & (d['interval'] == '10%')][[
#         'lower', 'upper']].values.mean(axis=1)
#     for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
#         width = _d[_d['quantity'] == 'width_mu']
#         frac = _d[_d['quantity'] == 'peri_drymass_frac']
#         ax.vlines(frac_10 * 100, 1/width['lower'], 1/width['upper'],
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+1, label='__nolegend__')
#         ax.hlines(1/width_10, frac['lower'] * 100, frac['upper'] * 100,
#                   lw=2, color=ppc_cmap_dict[_g], zorder=i+1, label='__nolegend__')

# k_range = [0.1]
# # k_range = np.linspace(0.05, 0.2, 5)
# phi_range = np.linspace(0.001, .15, 200)
# for _k in k_range:
#     theo_width = 4 * delta * (_k * (phi_range**-1 - 1) + 1) + 0.25
#     ax.plot(phi_range*100, 1/theo_width, lw=1,
#             label='model prediction, $k = 0.1$')

# ax.plot([], [], 'X', color='grey', label='mass spectrometry data')
# ax.plot([], [], 'P', lw=2, color=cor['primary_red'],
#         label='inferred from our measurements')
# ax.legend()
# ax.set_title(
#     r'predicted width assuming constant $k$ and $\ell/\omega \approx 4$')
# ax.set_ylim([0.25, 2.5])
# ax.set_xlim([0, 6])

# ax.set_xlabel('periplasmic protein drymass fraction [%]')
# ax.set_ylabel('width [µm]')

# plt.savefig('../../figures/mcmc/one-shot_revised_model_prediction.pdf',
#             bbox_inches='tight')
