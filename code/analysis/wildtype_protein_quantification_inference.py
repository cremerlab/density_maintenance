# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import cmdstanpy
import arviz
import corner
cor, pal = size.viz.matplotlib_style()

# Load data
calib_curve = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
protein_meas = pd.read_csv(
    '../../data/protein_quantification/wildtype_bradford_periplasmic_protein.csv')

# Load/compile inferential model
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/hierarchical_protein_quantification.stan')
# %%
# Assign indices
protein_meas['idx'] = protein_meas.groupby(['carbon_source']).ngroup() + 1
calib_curve['idx'] = calib_curve.groupby(
    ['replicate', 'protein_standard']).ngroup() + 1

# Generate the data dictionary
data_dict = {'J': protein_meas['idx'].max(),
             'K': calib_curve['idx'].max(),
             'N_meas': len(protein_meas),
             'N_calib': len(calib_curve),
             'meas_idx': protein_meas['idx'].values.astype(int),
             'calib_idx': calib_curve['idx'].values.astype(int),
             'concentration': calib_curve['protein_conc_ug_ml'].values.astype(float),
             'od_595nm_calib': calib_curve['od_595nm'].values.astype(float),
             'od_595nm_meas':  protein_meas['od_595nm'].values.astype(float),
             'od_600nm_meas': protein_meas['od_600nm'].values.astype(float),
             'dil_factor': protein_meas['dilution_factor'].values.astype(float),
             'ext_volume': protein_meas['extract_volume'].values.astype(float),
             'cult_volume': protein_meas['culture_volume'].values.astype(float),
             }

# Sample the model
samples = model.sample(data=data_dict, iter_warmup=6000, adapt_delta=0.99)
# %%
samples = arviz.from_cmdstanpy(samples)
# %%
fig = plt.figure(figsize=(8, 8))
_ = corner.corner(samples,
                  var_names=['slope',
                             'intercept',
                             'od_per_biomass_mu',
                             'calib_sigma'],
                  divergences=False,
                  smooth=1,
                  fig=fig,
                  divergences_kwargs={
                      'color': cor['light_blue'],
                      'ms': 4,
                      'markeredgewidth': 0})
fig.text(0, 1, 'wildtype protein quantification')
# plt.savefig('../../figures/wildtype_periplasmic_protein_diagnostics.pdf',
# bbox_inches='tight')
# %%
fig, ax = plt.subplots(1, 3, figsize=(7, 2))

# Plot the ppcs
calib_rep = samples.posterior.od_595nm_calib_rep.to_dataframe().reset_index()

for i in range(3):
    ind = np.where(calib_curve['idx'] == i+1)[0]
    concs = calib_curve[calib_curve['idx'] == i+1]
    ppc = calib_rep[calib_rep['od_595nm_calib_rep_dim_0'].isin(ind)]
    j = 0
    for g, d in ppc.groupby(['chain', 'draw']):
        if j % 10 == 0:
            ax[i].plot(concs['protein_conc_ug_ml'], d['od_595nm_calib_rep'],
                       '-', color=cor['primary_black'], lw=0.1)
        j += 1


i = 0
for g, d in calib_curve.groupby(['protein_standard', 'replicate']):
    ax[i].plot(d['protein_conc_ug_ml'], d['od_595nm'],
               '-o', color=cor['primary_red'])
    ax[i].set_title(f'{g[0]} standard, replicate {g[1]}')
    ax[i].set_xlabel('concentration [µg/mL]')
    ax[0].set_ylabel('OD$_{595nm}$')
    i += 1
