# %%
import numpy as np
import pandas as pd
import cmdstanpy
import arviz as az

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
bradford_prot_data['od_meas'] = bradford_prot_data['od_595nm'] * bradford_prot_data['dilution_factor'] * \
    bradford_prot_data['extraction_volume_ml'] / \
    (bradford_prot_data['od_600nm'] * bradford_prot_data['culture_volume_ml'])
bradford_prot_data['cond_idx'] = bradford_prot_data.groupby(
    ['carbon_source']).ngroup() + 1
# bradford_prot_data['brep_idx'] = bradford_prot_data.groupby(['date', 'carbon_source']).ngroup() + 1
bradford_cal_data['brep_idx'] = bradford_cal_data.groupby(
    ['replicate', 'protein_standard']).ngroup() + 1
growth_data['cond_idx'] = growth_data.groupby(['carbon_source']).ngroup() + 1
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
# , adapt_delta=0.99, max_treedepth=15)
_samples = model.sample(data=data_dict, max_treedepth=15, adapt_delta=0.95)  # , max_treedepth=15,
# adapt_delta=0.99, iter_sampling=2000, iter_warmup=1000)

# %%
samples = az.from_cmdstanpy(_samples)

#%%


# %%
cpb = samples.posterior.growth_cells_per_biomass.to_dataframe().reset_index()
cpb

# %%
rho = samples.posterior.peri_density.to_dataframe().reset_index()
rho

# %%
ppb = samples.posterior.prot_per_biomass.to_dataframe().reset_index()
ppb

#%%
cal = samples.posterior.cal_slope
cal
# %%
growth_mu_df = samples.posterior.growth_mu.to_dataframe().reset_index()
growth_mu_df.groupby(['growth_mu_dim_0'])['growth_mu'].mean()
growth_mu_df


#%%
width_mu = samples.posterior.width_mu.to_dataframe().reset_index()
width_mu