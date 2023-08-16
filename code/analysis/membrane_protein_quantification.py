# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()

# ##############################################################################
# DATASET LOADING
# ##############################################################################

# BCA calibration curve
bca_cal = pd.read_csv('../../data/bca_calibration_curve.csv')
bca_cal.sort_values(by='bsa_conc_ug_mL', inplace=True)
brad_cal = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
bca_meas = pd.read_csv('../../data/collated_BCA_measurements.csv')

# Membrane protein quantification
mem = bca[(bca['fraction'] == 'membrane') & (bca['overexpression'] == 'none')]
mem_agg = mem.groupby(['strain', 'overexpression', 'inducer_conc',
                      'carbon_source', 'biological_replicate']).mean().reset_index()
valid = []
for g, d in mem_agg.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc']):
    if len(d) >= 3:
        valid.append(d)
mem_valid = pd.concat(valid, sort=False)
mem_valid.dropna(inplace=True)
mem_valid['conv_factor'] = (mem_valid['od600nm'].values * mem_valid['culture_volume_mL'] /
                            (mem_valid['dilution_factor'].values * mem_valid['extraction_volume']))

# Periplasmic protein measurements
peri_data = pd.read_csv(
    '../../data/summaries/summarized_protein_measurements.csv')
peri_data = peri_data[(peri_data['strain'] == 'wildtype') & (peri_data['temperature_C'] == 37) & (
    peri_data['overexpression'] == 'none') & (peri_data['inducer_conc_ng_mL'] == 0)]
peri_data['conv_factor'] = peri_data['od_600nm'] * peri_data['dilution_factor'] * \
    peri_data['extraction_volume_mL'] / peri_data['culture_volume_mL']
# Growth rate data
growth = pd.read_csv('../../data/summaries/summarized_growth_measurements.csv')
growth = growth[(growth['strain'] == 'wildtype') & (growth['overexpression'] == 'none') &
                (growth['inducer_conc'] == 0) & (growth['temperature'] == 37)]

# Flow data
flow = pd.read_csv('../../data/summaries/flow_cytometry_counts.csv')
flow = flow[flow['strain'] == 'wildtype']


# Cell size data
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
size_data = size_data[(size_data['strain'] == 'wildtype') & (size_data['overexpression'] == 'none') & (
    size_data['inducer_conc'] == 0) & (size_data['temperature_C'] == 37) &
    (size_data['carbon_source'] != 'ezMOPS')]


# Prot per cell data
prot = pd.read_csv('../../data/literature/collated_total_protein.csv')


# ##############################################################################
# CONDITION MAPPING
# ##############################################################################
mapper = {'acetate': 1, 'sorbitol': 2, 'glycerol': 3,
          'glucose': 4, 'glucoseCAA': 5, 'LB': 6}

# For each dataset, apply the index labeling
mem_valid['idx'] = [mapper[s] for s in mem_valid['carbon_source'].values]
growth['idx'] = [mapper[s] for s in growth['carbon_source'].values]
flow['idx'] = [mapper[s] for s in flow['carbon_source'].values]
size_data['idx'] = [mapper[s] for s in size_data['carbon_source'].values]
peri_data['idx'] = [mapper[s] for s in peri_data['carbon_source'].values]

# %%
model = cmdstanpy.CmdStanModel(
    stan_file='stan_models/membrane_protein_quantification.stan')

# %%
data_dict = {
    #  Calibration curve
    'N_bca_cal': len(bca_cal),
    'N_brad_cal': len(brad_cal),
    'bca_std_conc': bca_cal['bsa_conc_ug_mL'].values,
    'brad_std_conc': brad_cal['protein_conc_ug_ml'].values,
    'bca_cal': bca_cal['od562nm'].values,
    'brad_cal': brad_cal['od_595nm'].values,

    # Membrane protein measurements
    'N_mem':  len(mem_valid),
    'J_mem': mem_valid['idx'].max(),
    'mem_idx': mem_valid['idx'].values,
    'mem_conv_factor': mem_valid['conv_factor'].values,
    'mem_od562nm': mem_valid['od562nm'].values,

    # Periplasmic protein measurement
    'N_peri':  len(peri_data),
    'J_peri': peri_data['idx'].max(),
    'peri_idx': peri_data['idx'].values,
    'peri_conv_factor': peri_data['conv_factor'].values,
    'peri_od595nm': peri_data['od_595nm'].values,

    # Flow
    'N_flow': len(flow),
    'J_flow': flow['idx'].max(),
    'flow_idx': flow['idx'].values,
    'cell_count': flow['cells_per_biomass'].values,

    # Growth
    'N_growth': len(growth),
    'J_growth': growth['idx'].max(),
    'growth_idx': growth['idx'].values,
    'growth_rates': growth['growth_rate_hr'].values,

    # Size data
    'N_size': len(size_data),
    'J_size': size_data['idx'].max(),
    'size_idx': size_data['idx'].values,
    'widths': size_data['width_median'].values,
    'lengths': size_data['length'].values,
    'volume': size_data['volume'].values,
    'aspect_ratio': size_data['aspect_ratio'].values,
    'surface_area': size_data['surface_area'].values,

    # Protein
    'N_prot': len(prot),
    'prot_growth_rate': prot['growth_rate_hr'].values,
    'prot_per_cell': prot['fg_protein_per_cell'].values
}

_samples = model.sample(data_dict)
samples = az.from_cmdstanpy(_samples)
# %%
cal_ppc = samples.posterior['bca_cal_ppc'].to_dataframe().reset_index()

fig, ax = plt.subplots()
for i, (g, d) in enumerate(cal_ppc.groupby(['chain', 'draw'])):
    if i % 2 == 0:
        ax.plot(cal['bsa_conc_ug_mL'].values,
                d['bca_cal_ppc'], 'k-', lw=0.1, alpha=0.15)

ax.plot(cal['bsa_conc_ug_mL'], cal['od562nm'], 'o', color=cor['primary_blue'])
ax.set_xscale('log')
ax.set_xlabel('BSA standard concentration [µg/mL]')
ax.set_ylabel('OD$_{562nm}$')

# %%
ppcs = samples.posterior['mem_ppc'].to_dataframe().reset_index()
fig, ax = plt.subplots(1, 1, figsize=(4, 6))
# ax.set_xlim(0, 5)
labs = []
for g, d in mem_valid.groupby(['idx', 'carbon_source', 'overexpression', 'inducer_conc']):
    labs.append(f'{g[1:]}')
    inds = np.where(mem_valid['idx'].values == g[0])[0]
    _ppcs = ppcs[ppcs['mem_ppc_dim_0'].isin(inds)]
    ax.plot(d['od562nm'],  g[0] * np.ones(len(d)),
            'o', color=cor['primary_red'], zorder=1000)
    ax.plot(_ppcs['mem_ppc'], g[0] * np.ones(len(_ppcs)), '|',
            markeredgecolor='k', markeredgewidth=0.1, alpha=0.1)


ax.set_yticks(mem_valid['idx'].unique())
ax.set_yticklabels(labs)


# ax.set_yticklabels(mem_valid.groupby(['idx']))
# %%
ms = pd.read_csv('../../data/literature/collated_mass_fractions_empirics.csv')
phi_mem_post = samples.posterior['phi_mem'].to_dataframe().reset_index()
lam_mu_post = samples.posterior['lam_mu'].to_dataframe().reset_index()
fig, ax = plt.subplots(1, 1, figsize=(3, 2))

for g, d in ms[ms['localization'] == 'membrane'].groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['mass_frac'], **fmt)

for k, v in mapper.items():
    phi_mem = phi_mem_post[phi_mem_post['phi_mem_dim_0'] == v - 1]
    lam_mu = lam_mu_post[lam_mu_post['lam_mu_dim_0'] == v - 1]
    med_phi = np.median(phi_mem['phi_mem'].values)
    med_lam = np.median(lam_mu['lam_mu'].values)
    phi_percs = np.percentile(phi_mem['phi_mem'].values, [97.5, 2.5])
    lam_percs = np.percentile(lam_mu['lam_mu'].values, [97.5, 2.5])

    ax.vlines(med_lam, phi_percs[0], phi_percs[1],
              linewidth=1, color=cor['primary_blue'])
    ax.hlines(med_phi, lam_percs[0], lam_percs[1],
              linewidth=1, color=cor['primary_blue'])
    ax.plot(med_lam, med_phi,  'o', markerfacecolor='w', markeredgecolor=cor['primary_blue'],
            markeredgewidth=0.5, ms=4.5)

ax.set_ylim([0, 0.25])
ax.set_ylabel('$\phi_{mem}$', fontsize=8)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)

# %%
rho_mem_post = samples.posterior['rho_mem'].to_dataframe().reset_index()
lam_mu_post = samples.posterior['lam_mu'].to_dataframe().reset_index()
fig, ax = plt.subplots(1, 1, figsize=(3, 2))

for g, d in ms[ms['localization'] == 'membrane'].groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['mass_fg'] / (2 * d['surface_area']), **fmt)

for k, v in mapper.items():
    rho_mem = rho_mem_post[rho_mem_post['rho_mem_dim_0'] == v - 1]
    lam_mu = lam_mu_post[lam_mu_post['lam_mu_dim_0'] == v - 1]
    med_rho = np.median(rho_mem['rho_mem'].values)
    med_lam = np.median(lam_mu['lam_mu'].values)
    rho_percs = np.percentile(rho_mem['rho_mem'].values, [97.5, 2.5])
    lam_percs = np.percentile(lam_mu['lam_mu'].values, [97.5, 2.5])

    ax.vlines(med_lam, rho_percs[0], phi_percs[1],
              linewidth=1, color=cor['primary_blue'])
    ax.hlines(med_rho, lam_percs[0], lam_percs[1],
              linewidth=1, color=cor['primary_blue'])
    ax.plot(med_lam, med_rho,  'o', markerfacecolor='w', markeredgecolor=cor['primary_blue'],
            markeredgewidth=1, ms=4.5)

# ax.set_ylim([0, 0.25])
ax.set_ylabel('$\phi_{mem}$', fontsize=8)
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=8)
