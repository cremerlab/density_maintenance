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
bradford_prot_data['od_meas'] = (bradford_prot_data['od_595nm'] * bradford_prot_data['extraction_volume_ml'] *
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
size_data['cond_idx'] = size_data.groupby(['carbon_source']).ngroup() + 1

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
    'prot_od': bradford_prot_data['od_meas'].values.astype(float),
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
        if type(g) != list:
            g = [g]
        # Compute the percentiles for different quantities
        for q in quantity:
            lower = np.percentile(d[f'{q}'].values, lower_bounds)
            upper = np.percentile(d[f'{q}'].values, upper_bounds)
            _df = pd.DataFrame(np.array([lower, upper]).T, columns=[
                               'lower', 'upper'])
            _df['quantity'] = q
            _df['interval'] = interval_labels

            # Add the grouping informaton
            for i, _g in enumerate(g):
                _df[groupby[i]] = _g
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
# %%
