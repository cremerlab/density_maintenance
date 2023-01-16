# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
import size.viz
import seaborn as sns
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

# Restrict the datasets to wildtype
bradford_prot_data = bradford_prot_data[
    (bradford_prot_data['strain'] == 'wildtype') &
    (bradford_prot_data['inducer_conc_ng_ml'] == 0)
]
growth_data = growth_data[growth_data['strain'] == 'wildtype']
flow_data = flow_data[flow_data['strain'] == 'wildtype']
flow_data = flow_data.groupby(['carbon_source', 'date'])[
    'cells_per_biomass'].mean().reset_index()
size_data = size_data[(size_data['strain'] == 'wildtype')
                      & (size_data['inducer_conc'] == 0) &
                      (size_data['temperature_C'] == 37) &
                      (size_data['carbon_source'] != 'ezMOPS')]

# Compile the inferential model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/one-shot_semicomplete_inference.stan')
# %%
# Compute the necessary properties and idx groups
bradford_prot_data['conv_factor'] = (bradford_prot_data['extraction_volume_ml'] *
                                 bradford_prot_data['dilution_factor']) / (bradford_prot_data['od_600nm'] * bradford_prot_data['culture_volume_ml'])

bradford_prot_data['cond_idx'] = bradford_prot_data.groupby(
    ['carbon_source']).ngroup() + 1
bradford_cal_data['brep_idx'] = bradford_cal_data.groupby(
    ['replicate', 'protein_standard']).ngroup() + 1

# Hardcode the indexing for the growth data conditions
growth_data['cond_idx'] = growth_data.groupby(
    ['carbon_source'], sort=False).ngroup()

growth_data.loc[growth_data['carbon_source'] == 'LB', 'cond_idx'] = 6


growth_data['brep_idx'] = growth_data.groupby(
    ['carbon_source', 'date', 'run_no']).ngroup() + 1


# Map flow cytometry data to appropriate carbon source index
carb_idx_dict = {g[0]: g[1] for g, _ in growth_data[growth_data['elapsed_time_hr'] == 0].groupby([
    'carbon_source', 'cond_idx'])}
flow_data['growth_idx'] = [carb_idx_dict[i]
                           for i in flow_data['carbon_source'].values]
carb_idx = [v for _, v in carb_idx_dict.items()]

# Add identifying information to size data
size_data['cond_idx'] = size_data.groupby(['carbon_source']).ngroup()
size_data.loc[size_data['carbon_source'] == 'LB', 'cond_idx'] = 6

# Define the data dictionary
data_dict = {
    # Bradford calibration curve inputs
    'J_cal_brep': bradford_cal_data['brep_idx'].max(),
    'N_cal_meas': len(bradford_cal_data),
    'cal_idx': bradford_cal_data['brep_idx'].values.astype(int),
    'concentration': bradford_cal_data['protein_conc_ug_ml'].values.astype(float),
    'cal_od': bradford_cal_data['od_595nm'].values.astype(float),

    # Protein quantification inputs
    'J_prot_cond': bradford_prot_data['cond_idx'].max(),
    'N_prot_meas': len(bradford_prot_data),
    'prot_idx': bradford_prot_data['cond_idx'].values.astype(int),
    'prot_od': bradford_prot_data['od_595nm'].values.astype(float),
    'od_conv_factor': bradford_prot_data['conv_factor'].values.astype(float),
    'mean_prot_od': bradford_prot_data.groupby(['cond_idx'])['od_meas'].mean().astype(float),
    'std_prot_od': bradford_prot_data.groupby(['cond_idx'])['od_meas'].std().astype(float),

    # Growth curve inputs
    'J_growth_cond': growth_data['cond_idx'].max(),
    'J_growth_brep': growth_data['brep_idx'].max(),
    'N_growth_meas': len(growth_data),
    'growth_cond_idx': growth_data.groupby(['brep_idx'])['cond_idx'].min().astype(int),
    'growth_brep_idx': growth_data['brep_idx'].values.astype(int),
    'growth_time': growth_data['elapsed_time_hr'].values.astype(float),
    'growth_od': growth_data['od_600nm'].values.astype(float),

    # Flow cytometry inputs
    'N_cell_meas': len(flow_data),
    'cells_per_biomass': flow_data['cells_per_biomass'].values.astype(float),

    # Size inputs
    'J_size_cond': size_data['cond_idx'].max(),
    'N_size_meas': len(size_data),
    'size_idx': size_data['cond_idx'].values.astype(int),
    'mean_width': size_data['width_median'].values.astype(float),
    'mean_length': size_data['length'].values.astype(float),
    'mean_vol': size_data['volume'].values.astype(float),
    'mean_peri_vol': size_data['periplasm_volume'].values.astype(float),
    'mean_peri_vol_frac': size_data['periplasm_volume_fraction'].values.astype(float),
    'mean_sa': size_data['surface_area'].values.astype(float),
    'mean_sav': size_data['surface_to_volume'].values.astype(float),
    'mean_width_mean': size_data.groupby(['cond_idx'])['width_median'].mean().astype(float),
    'mean_width_std': size_data.groupby(['cond_idx'])['width_median'].std().astype(float),
    'mean_length_mean': size_data.groupby(['cond_idx'])['length'].mean().astype(float),
    'mean_length_std': size_data.groupby(['cond_idx'])['length'].std().astype(float),
    'mean_vol_mean': size_data.groupby(['cond_idx'])['volume'].mean().astype(float),
    'mean_vol_std': size_data.groupby(['cond_idx'])['volume'].std().astype(float),
    'mean_peri_vol_mean': size_data.groupby(['cond_idx'])['periplasm_volume'].mean().astype(float),
    'mean_peri_vol_std': size_data.groupby(['cond_idx'])['periplasm_volume'].std().astype(float),
    'mean_peri_vol_frac_mean': size_data.groupby(['cond_idx'])['periplasm_volume_fraction'].mean().astype(float),
    'mean_peri_vol_frac_std': size_data.groupby(['cond_idx'])['periplasm_volume_fraction'].std().astype(float),
    'mean_sa_mean': size_data.groupby(['cond_idx'])['surface_area'].mean().astype(float),
    'mean_sa_std': size_data.groupby(['cond_idx'])['surface_area'].std().astype(float),
    'mean_sav_mean': size_data.groupby(['cond_idx'])['surface_to_volume'].mean().astype(float),
    'mean_sav_std': size_data.groupby(['cond_idx'])['surface_to_volume'].std().astype(float),

    # Linking indices
    'cells_per_biomass_growth_idx': flow_data['growth_idx'].values.astype(int),
    'prot_cond_map': np.array([ind for ind in carb_idx if ind in bradford_prot_data['cond_idx'].values]).astype(int),

}

# %%
# Sample the model
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
# Calibration curve ppc

# Unpack the calbiration curve ppc
concs = bradford_cal_data['protein_conc_ug_ml'].values
cal_ppc_df = samples.posterior.od595_calib_rep.to_dataframe().reset_index()
for i, c in enumerate(concs):
    cal_ppc_df.loc[cal_ppc_df['od595_calib_rep_dim_0']
                   == i, 'protein_conc_ug_ml'] = c


def compute_percentiles(df,
                        quantity,
                        groupby,
                        lower_bounds=[5, 10, 15, 20, 25, 30, 35, 40, 45, ],
                        # [0.5, 2.5, 12.5, 25, 37.5, 45, 49.5],
                        upper_bounds=[95, 90, 85, 80, 75, 70, 65, 60, 55],
                        # [99.5, 97.5, 87.5, 75, 62.5, 55, 50.5],
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


# %%
# PPC for calibration curve
cal_ppc_percs = compute_percentiles(
    cal_ppc_df, 'od595_calib_rep', 'protein_conc_ug_ml')
cal_ppc_percs

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
ppc_cmap = sns.color_palette('Greys_r', n_colors=len(
    cal_ppc_percs['interval'].unique()) + 2).as_hex()
ppc_cmap_dict = {i: c for i, c in zip(
    [f'{k}%' for k in np.arange(10, 100, 10)], ppc_cmap)}
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
plt.savefig('../../figures/mcmc/one-shot_calibration_curve_ppc.pdf')
# %%
# PPC for bradford assay
prot_ppc_df = samples.posterior.od595_per_biomass_rep.to_dataframe().reset_index()
# Map the index to the carbon source
_prot_data = bradford_prot_data[['carbon_source', 'cond_idx']]
_prot_data['n'] = np.arange(len(_prot_data))
for k, v in zip(_prot_data['carbon_source'].values, _prot_data['n'].values):
    prot_ppc_df.loc[prot_ppc_df['od595_per_biomass_rep_dim_0']
                    == v, 'carbon_source'] = k

prot_ppc_df

# Compute the percentiles
prot_ppc_percs = compute_percentiles(
    prot_ppc_df, 'od595_per_biomass_rep', 'carbon_source')

# Plot the PPC under the measurements
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
carb_loc = {'acetate': 1, 'sorbitol': 2, 'glycerol': 3,
            'glucose': 4, 'glucoseCAA': 5, 'LB': 6}
for i, (g, d) in enumerate(prot_ppc_percs.groupby(['carbon_source', 'interval'], sort=False)):
    ax.hlines(carb_loc[g[0]], d['lower'], d['upper'],
              color=ppc_cmap_dict[g[1]], lw=15, zorder=i)

for g, d in bradford_prot_data.groupby(['carbon_source']):
    ax.plot(d['od_meas'], carb_loc[g] + np.random.normal(0, 0.01,
            len(d)), 'o', ms=3, color=cor['primary_red'], zorder=1000)

_ = ax.set_yticks(list(carb_loc.values())[:-1])
_ = ax.set_yticklabels(list(carb_loc.keys())[:-1])
_ = ax.set_xlabel('OD$_{595nm}$ per biomass')
ax.set_title('Periplasmic protein quantification measurements')
plt.tight_layout()
plt.savefig('../../figures/mcmc/one-shot_periprot_quantification_ppc.pdf')
# %%
# %%
cpb_ppc_df = samples.posterior[[
    'cells_per_biomass_rep']].to_dataframe().reset_index()

# Map conditions to carbon_source
flow_data['n'] = np.arange(len(flow_data))
for dim, carb in zip(flow_data['n'].values, flow_data['carbon_source'].values):
    cpb_ppc_df.loc[cpb_ppc_df['cells_per_biomass_rep_dim_0']
                   == dim, 'carbon_source'] = carb

# Compute the percentiles
cpb_perc_df = compute_percentiles(
    cpb_ppc_df, 'cells_per_biomass_rep', 'carbon_source')
fig, ax = plt.subplots(1, 1, figsize=(4, 2))

for i, (g, d) in enumerate(cpb_perc_df.groupby(['carbon_source', 'interval'], sort=False)):
    ax.hlines(carb_loc[g[0]], d['lower'], d['upper'],
              zorder=i, lw=15, color=ppc_cmap_dict[g[1]])
ax.set_xscale('log')

for g, d in flow_data.groupby(['carbon_source']):
    ax.plot(d['cells_per_biomass'], carb_loc[g] + np.random.normal(0, 0.1, len(d)), 'o',
            color=cor['primary_red'], ms=4, zorder=1000)
ax.set_xlim([1E8, 5E9])
ax.set_yticks([2, 3, 4, 5, 6])
ax.set_yticklabels(['sorbitol', 'glycerol', 'glucose', 'glucoseCAA', 'LB'])
ax.set_xlabel('cells per biomass')
ax.set_title('Flow cytometry event counts')
plt.tight_layout()
plt.savefig('../../figures/mcmc/one-shot_flow_cytometry_ppc.pdf')
# %%
# Plot the fit to the flow data
cpb_fit_df = samples.posterior[[
    'growth_cells_per_biomass']].to_dataframe().reset_index()
for carb, dim in carb_idx_dict.items():
    cpb_fit_df.loc[cpb_fit_df['growth_cells_per_biomass_dim_0']
                   == dim-1, 'carbon_source'] = carb

_df = samples.posterior[[
    'growth_mu']].to_dataframe().reset_index()
cpb_fit_df['growth_mu'] = _df['growth_mu'].values

# Compute the percentiles
cpb_fit_percs = compute_percentiles(
    cpb_fit_df, ['growth_cells_per_biomass', 'growth_mu'], 'carbon_source')
cpb_fit_percs

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
for g, d in cpb_fit_percs.groupby(['carbon_source'], sort=False):
    lam = d[d['quantity'] == 'growth_mu']
    cpb = d[d['quantity'] == 'growth_cells_per_biomass']
    perc_10_cpb = cpb[cpb['interval'] == '10%']
    perc_10_lam = lam[lam['interval'] == '10%']
    x = np.mean([perc_10_lam['lower'].values[0],
                perc_10_lam['upper'].values[0]])
    y = np.mean([perc_10_cpb['lower'].values[0],
                perc_10_cpb['upper'].values[0]])
    for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
        _lam = _d[_d['quantity'] == 'growth_mu']
        _cpb = _d[_d['quantity'] == 'growth_cells_per_biomass']
        ax.vlines(x, _cpb['lower'], _cpb['upper'], lw=3,
                  zorder=i+1, color=ppc_cmap_dict[_g])
        ax.hlines(y, _lam['lower'], _lam['upper'], lw=3,
                  zorder=i+1, color=ppc_cmap_dict[_g])

    if g in flow_data['carbon_source'].values:
        _data = flow_data[flow_data['carbon_source'] == g]
        ax.plot(x * np.ones(len(_data)),
                _data['cells_per_biomass'], 'o', ms=5, color=cor['primary_red'], zorder=1000)

# compute and plot the main fit
fit_params = samples.posterior[[
    'k_cells_per_biomass_tilde', 'beta_0_tilde']].to_dataframe().reset_index()
lam_range = np.linspace(0.1, 2.5, 200)
fit = 1E9 * np.exp(fit_params['beta_0_tilde'].mean() -
                   fit_params['k_cells_per_biomass_tilde'].mean() * lam_range)
ax.plot(lam_range, fit, '-', lw=1, color=cor['primary_blue'])
ax.set_xlabel('growth rate [hr$^{-1}$]')
ax.set_ylabel('cells per biomass')
ax.set_title('Exponential $\lambda$-dependence on cell density')
plt.tight_layout()
plt.savefig('../../figures/mcmc/one-shot_cell_density_lam_dependence.pdf')

# %%
# Plot the growth ppcs
lam_ppc_df = samples.posterior.growth_od_rep.to_dataframe().reset_index()

# Map time and carbon sources to everything
carbs = growth_data['carbon_source'].values
breps = growth_data['brep_idx'].values
elapsed_time = growth_data['elapsed_time_hr'].values
dims = np.arange(len(growth_data))
for dim, c, b, t in zip(dims, carbs, breps, elapsed_time):
    lam_ppc_df.loc[lam_ppc_df['growth_od_rep_dim_0']
                   == dim, 'carbon_source'] = c
    lam_ppc_df.loc[lam_ppc_df['growth_od_rep_dim_0'] == dim, 'brep'] = int(b)
    lam_ppc_df.loc[lam_ppc_df['growth_od_rep_dim_0'] == dim, 'time'] = t
lam_ppc_percs = compute_percentiles(lam_ppc_df, 'growth_od_rep', [
                                    'carbon_source', 'brep', 'time'])
# %%
for g, d in growth_data.groupby(['carbon_source']):
    n_reps = len(d[d['elapsed_time_hr'] == 0])
    cols = 3
    rows = int(np.ceil(n_reps / cols))
    ex = cols * rows - n_reps
    fig, ax = plt.subplots(rows, cols, figsize=(8, 2 * rows))

    # Turn off unpopulated axes
    if ex != 0:
        for a in ax.ravel()[-ex:]:
            a.axis(False)
        ax = ax.ravel()[:-ex]
    else:
        ax = ax.ravel()
    for a in ax:
        a.set_xlabel('elapsed time [hr]')
        a.set_ylabel('log OD$_{600nm}$')

    # get the ppc
    _ppc = lam_ppc_percs[(lam_ppc_percs['carbon_source'] == g)]
    ind = 0
    for i, (_g, _d) in enumerate(_ppc.groupby(['brep', 'interval'], sort=False)):
        if i == 0:
            _brep = _g[0]
        if _g[0] != _brep:
            ind += 1
            _brep = _g[0]

        ax[ind].fill_between(_d['time'], np.log(_d['lower']), np.log(_d['upper']),
                             color=ppc_cmap_dict[_g[1]], zorder=i+1)

    for i, (_g, _d) in enumerate(d.groupby(['brep_idx'])):
        ax[i].plot(_d['elapsed_time_hr'], np.log(_d['od_600nm']),
                   '-o', color=cor['primary_red'], ms=4, zorder=1000, lw=1)
    plt.suptitle(f'wildtype {g} growth curves')
    plt.tight_layout()
    plt.savefig(f'../../figures/mcmc/one-shot_{g}_growth_curves_ppc.pdf')
    plt.close()

# %%
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

# %%


samples.posterior.peri_density.to_dataframe().reset_index().groupby([
    'peri_density_dim_0']).median()
