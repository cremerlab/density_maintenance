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

# Restrict the datasets to wildtype
bradford_prot_data = bradford_prot_data[
    (bradford_prot_data['strain'] == 'wildtype') &
    (bradford_prot_data['inducer_conc_ng_ml'] == 0)
]
growth_data = growth_data[growth_data['strain'] == 'wildtype']

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
}


# %%
# Sample the model
_samples = model.sample(data=data_dict, adapt_delta=0.95, max_treedepth=15)

# %%
samples = az.from_cmdstanpy(samples)
# %%
growth_mu_df = samples.posterior.growth_mu.to_dataframe().reset_index()
growth_mu_df.groupby(['growth_mu_dim_0'])['growth_mu'].mean()


growth_data
