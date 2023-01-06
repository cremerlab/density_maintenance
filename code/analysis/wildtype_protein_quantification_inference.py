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
data_dict = {'J_cond': protein_meas['idx'].max(),
             'J_calib': calib_curve['idx'].max(),
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

# %% save hyperparameter distributions
prot_per_biomass = samples.posterior['prot_per_biomass'].to_dataframe(
).reset_index()
for g, d in protein_meas.groupby(['carbon_source', 'idx']):
    prot_per_biomass.loc[prot_per_biomass['prot_per_biomass_dim_0']
                         == (g[1] - 1), 'carbon_source'] = g[0]
prot_per_biomass['strain'] = 'wildtype'
prot_per_biomass = prot_per_biomass[[
    'carbon_source', 'prot_per_biomass', 'strain']]
prot_per_biomass.to_csv(
    '../../data/protein_quantification/mcmc/wildtype_protein_per_biomass_hyperparameter_samples.csv', index=False)


# %% Compute summary statistics
percs = [2.5, 12.5, 25, 45, 55, 75, 87.5, 97.5]
perc_cols = ['2.5%', '12.5%', '25%', '45%', '55%', '75%', '87.5%', '97.5%']
summary_df = pd.DataFrame([])
for g, d in prot_per_biomass.groupby(['carbon_source']):
    _percs = pd.DataFrame([np.percentile(d['prot_per_biomass'], percs)],
                          columns=perc_cols)
    _percs['mean'] = d['prot_per_biomass'].mean()
    _percs['median'] = d['prot_per_biomass'].median()
    _percs['carbon_source'] = g
    _percs['strain'] = 'wildtype'
    _percs['parameter'] = 'ug protein per OD mL'
    summary_df = pd.concat([summary_df, _percs])
summary_df.to_csv(
    '../../data/protein_quantification/mcmc/wildtype_protein_per_biomass_hyperparameter_summary.csv', index=False)

# %%
fig = plt.figure(figsize=(8, 8))
_ = corner.corner(samples,
                  var_names=['slope',
                             'intercept',
                             'od_per_biomass_mu',
                             ],
                  divergences=False,
                  smooth=1,
                  fig=fig,
                  divergences_kwargs={
                      'color': cor['light_blue'],
                      'ms': 4,
                      'markeredgewidth': 0})
fig.text(0, 1, 'wildtype protein quantification')
plt.savefig('../../figures/mcmc/protein_diagnostics/wildtype_periplasmic_protein_diagnostics.pdf',
            bbox_inches='tight')
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
plt.savefig('../../figures/mcmc/protein_diagnostics/bradford_calibration_ppc.pdf',
            bbox_inches='tight')


# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 2))
ax.set_yticks(protein_meas['idx'].unique())
ax.set_yticklabels(protein_meas['carbon_source'].unique())

# Plot ppc for OD per biomass
od595_rep = samples.posterior.od_per_biomass_rep.to_dataframe().reset_index()
for i in protein_meas['idx'].unique():
    ind = np.where(protein_meas['idx'] == i)[0]
    ppc = od595_rep[od595_rep['od_per_biomass_rep_dim_0'].isin(ind)]
    for g, d in ppc.groupby(['chain', 'draw']):
        if j % 10 == 0:
            ax.plot(d['od_per_biomass_rep'], np.ones(len(d)) * i + np.random.normal(0, 0.05, len(d)),
                    '.', color=cor['primary_black'], lw=0.1, markeredgewidth=0, alpha=0.5, ms=2)
for g, d in protein_meas.groupby(['idx']):
    d['resc'] = d['od_595nm'] * d['dilution_factor'] * \
        (d['extract_volume'].values) / \
        (d['culture_volume'] * d['od_600nm'])
    ax.plot(d['resc'], np.ones(len(d)) * g + np.random.normal(0, 0.05, len(d)), 'o', ms=3,
            color=cor['primary_red'], markeredgewidth=0.5)
ax.set_xlabel('OD$_{595nm}$ per OD$_{600nm}$ mL', fontsize=6)
plt.savefig('../../figures/mcmc/protein_diagnostics/wildtype_OD_per_biomass_ppc.pdf',
            bbox_inches='tight')
