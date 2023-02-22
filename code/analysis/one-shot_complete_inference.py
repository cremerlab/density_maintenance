# %%
import numpy as np
import pandas as pd
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
    stan_file='./stan_models/one-shot_complete_inference.stan')

# %%
# Load the various datasets
cal_data = pd.read_csv(
    '../../data/protein_quantification/bradford_calibration_curve.csv')
brad_data = pd.read_csv(
    '../../data/protein_quantification/bradford_periplasmic_protein_v2.csv')
brad_data['od_600nm_true'] = brad_data['od_600nm']
# brad_data['od_600nm'] *= 0.78
brad_data.loc[brad_data['od_600nm'] >= 0.45, 'od_600nm_true'] = np.exp(
    1.26 * np.log(brad_data[brad_data['od_600nm'] >= 0.45]['od_600nm'].values) + 0.25)
brad_data = brad_data[brad_data['od_600nm'] <= 0.45]
brad_data = brad_data[~((brad_data['overexpression'] != 'none') & (
    brad_data['inducer_conc_ng_mL'] == 0))]
size_data = pd.read_csv(
    '../../data/summaries/summarized_size_measurements.csv')
biomass_data = pd.read_csv(
    '../../data/literature/Basan2015/Basan2015_drymass_protein_cellcount.csv')
# %%
# Keep only the bradford data with more than two replicates
brad_data = pd.concat([d for _, d in brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']) if len(d) > 2], sort=False)
brad_data = brad_data[brad_data['strain'].isin(
    ['wildtype',  'malE-rbsB-fliC-KO'])]
size_data = size_data[size_data['temperature_C'] == 37]
size_data = size_data[size_data['strain'].isin(
    ['wildtype', 'malE-rbsB-fliC-KO'])]
size_data = pd.concat([d for _, d in size_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc']) if len(d) > 2], sort=False)
size_data = size_data[~((size_data['overexpression'] != 'none') & (
    size_data['inducer_conc'] == 0))]

# %%
# Filter, label, and transform bradford data
brad_data = brad_data[brad_data['strain'].isin(
    ['wildtype',  'malE-rbsB-fliC-KO'])]
brad_data = brad_data[~((brad_data['overexpression'] != 'none') & (
    brad_data['inducer_conc_ng_mL'] == 0))]
brad_data['cond_idx'] = brad_data.groupby(
    ['strain', 'carbon_source', 'overexpression', 'inducer_conc_ng_mL']).ngroup() + 1
brad_data['conv_factor'] = brad_data['dilution_factor'] * \
    brad_data['extraction_volume_mL'] / \
    (brad_data['culture_volume_mL'])
brad_data['od_per_biomass'] = brad_data['od_595nm']

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
             'brad_od595': brad_data['od_595nm'].values.astype(float),
             'brad_od600': brad_data['od_600nm_true'].values.astype(float),
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
_samples = model.sample(data_dict, adapt_delta=0.95, iter_sampling=5000)
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
ax.set_xlabel('OD$_{595nm}$')

# %%
# Biomass ppc
biomass_ppc = samples.posterior.biomass_rep.to_dataframe()
biomass_ppc['idx'] = 1
biomass_percs = size.viz.compute_percentiles(biomass_ppc, 'biomass_rep', 'idx')

fig, ax = plt.subplots(1, 1, figsize=(4, 1))

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

# width_ppc.dropna(inplace=True)
width_perc = size.viz.compute_percentiles(width_ppc, 'width_mu', [
                                          'link_idx', 'strain', 'carbon_source', 'overexpression', 'inducer_conc'])


# %%
# Mass spec data
mass_fracs = pd.read_csv('../../data/literature/compiled_mass_fractions.csv')
mass_fracs = mass_fracs[mass_fracs['periplasm']]
dry_frac = 0.3
prot_frac = 0.55
density = 1.1

# %%
# Do the proper classification
# genes = pd.read_csv('../../data/literature/genes_classification_all.csv')
# _genes = genes[genes['location'].isin(['IM', 'OM', 'PE', 'LPO'])]
# mass_fracs = mass_fracs[mass_fracs['gene_name'].isin(_genes['gene'].unique())]
growth = pd.read_csv('../../data/summaries/summarized_growth_measurements.csv')
growth = growth[(growth['strain'] == 'wildtype') & (
    growth['overexpression'] == 'none') & (growth['inducer_conc_ng_ml'] == 0)]
growth = growth.groupby(['carbon_source']).mean().reset_index()
wt_size = size_data[(size_data['strain'] == 'wildtype') & (
    size_data['overexpression'] == 'none') & (size_data['inducer_conc'] == 0)]
for g, d in growth.groupby(['carbon_source', 'growth_rate_hr']):
    wt_size.loc[wt_size['carbon_source'] == g[0], 'growth_mu'] = g[1]

# Determine the simple relations
w_popt = scipy.stats.linregress(
    wt_size['growth_mu'], wt_size['width_median'])
ell_popt = scipy.stats.linregress(
    wt_size['growth_mu'], np.log(wt_size['length']))
peri_vol_popt = scipy.stats.linregress(
    wt_size['growth_mu'], np.log(wt_size['periplasm_volume']))
vol_popt = scipy.stats.linregress(
    wt_size['growth_mu'], np.log(wt_size['volume']))
sav_popt = scipy.stats.linregress(
    wt_size['growth_mu'], np.log(wt_size['surface_to_volume']))

# Compute the periplasmic protein density
mass_fracs['width'] = w_popt[0] * \
    mass_fracs['growth_rate_hr'] + w_popt[1]
mass_fracs['length'] = np.exp(
    ell_popt[0] * mass_fracs['growth_rate_hr'] + ell_popt[1])
mass_fracs['peri_vol'] = np.exp(
    peri_vol_popt[0] * mass_fracs['growth_rate_hr'] + peri_vol_popt[1])
mass_fracs['volume'] = np.exp(
    vol_popt[0] * mass_fracs['growth_rate_hr'] + vol_popt[1])
mass_fracs['sav'] = np.exp(
    sav_popt[0] * mass_fracs['growth_rate_hr'] + sav_popt[1])
# mass_fracs['peri_volume'] = size.analytical.surface_area(mass_fracs['length'], mass_fracs['width']) * 0.025
mass_fracs['tot_protein'] = density * \
    dry_frac * prot_frac * mass_fracs['volume']
mass_fracs['peri_protein'] = mass_fracs['mass_frac'] * \
    mass_fracs['tot_protein']
mass_fracs['rho_peri'] = (
    mass_fracs['peri_protein'] * 1E3) / mass_fracs['peri_vol']
mass_fracs['biomass_frac'] = mass_fracs['peri_protein'] / \
    (density * dry_frac * mass_fracs['volume'])
mass_fracs = mass_fracs.groupby(
    ['dataset_name', 'condition', 'growth_rate_hr', 'width']).sum().reset_index()


# %%
# Look only at the no inducer cases
noind_width = width_perc[(width_perc['inducer_conc'] == 0) & (width_perc['overexpression'] == 'none') & (
    width_perc['strain'].isin(['wildtype', 'malE-rbsB-KO', 'malE-rbsB-fliC-KO']))]
noind_frac = frac_perc[(frac_perc['inducer_conc'] == 0) & (frac_perc['overexpression'] == 'none') & (
    frac_perc['strain'].isin(['wildtype', 'malE-rbsB-KO', 'malE-rbsB-fliC-KO']))]

ko_ppc_cmap = {k: c for k, c in zip(frac_perc['interval'].unique(
), sns.color_palette('Blues_r', n_colors=3 + len(frac_perc['interval'].unique())))}
wt_ppc_cmap = {k: c for k, c in zip(frac_perc['interval'].unique(
), sns.color_palette('Greens_r', n_colors=3 + len(frac_perc['interval'].unique())))}
rbs_ppc_cmap = {k: c for k, c in zip(frac_perc['interval'].unique(
), sns.color_palette('Oranges_r', n_colors=3 + len(frac_perc['interval'].unique())))}

mal_ppc_cmap = {k: c for k, c in zip(frac_perc['interval'].unique(
), sns.color_palette('Purples_r', n_colors=3 + len(frac_perc['interval'].unique())))}
lac_ppc_cmap = {k: c for k, c in zip(frac_perc['interval'].unique(
), sns.color_palette('Greys_r', n_colors=3 + len(frac_perc['interval'].unique())))}

cmaps = {'wildtype': wt_ppc_cmap, 'malE-rbsB-fliC-KO': ko_ppc_cmap,
         'rbsB': rbs_ppc_cmap, 'malE': mal_ppc_cmap,
         'lacZ': lac_ppc_cmap}

width_inc = 0.3
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
for g, d in noind_frac.groupby(['strain', 'carbon_source']):
    _noind_width = noind_width[(noind_width['strain'] == g[0]) & (
        noind_width['carbon_source'] == g[1])]

    # find the x,y positions
    width10 = _noind_width[_noind_width['interval']
                           == '10%'][['lower', 'upper']].values.mean()
    frac10 = d[d['interval'] == '10%'][['lower', 'upper']].values.mean()
    ax.plot(frac10, width10, 'o', markeredgecolor='white',
            ms=4, markeredgewidth=0.5, color=cmaps[g[0]]['10%'], zorder=1000)
    for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
        ax.hlines(width10, _d['lower'], _d['upper'], color=cmaps[g[0]][_g],
                  zorder=i+1, label='__nolegend__', lw=0.1 + i * width_inc)

    for i, (_g, _d) in enumerate(_noind_width.groupby(['interval'], sort=False)):
        ax.vlines(frac10, _d['lower'], _d['upper'], color=cmaps[g[0]][_g],
                  zorder=i + 1, label='__nolegend__', lw=0.1 + i * width_inc)

ind_width = width_perc[(width_perc['overexpression'] != 'none') & (width_perc['inducer_conc'] > 0) & (
    width_perc['strain'].isin(['wildtype', 'malE-rbsB-fliC-KO']))]
ind_frac = frac_perc[(frac_perc['overexpression'] != 'none') & (frac_perc['inducer_conc'] > 0) & (
    frac_perc['strain'].isin(['wildtype', 'malE-rbsB-fliC-KO']))]


for g, d in ind_frac.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc']):
    # if g[0] == 'wildtype':
    # continue
    # if (g[2] != 'rbsB'):
    _ind_width = ind_width[(ind_width['strain'] == g[0]) & (
        ind_width['carbon_source'] == g[1]) & (ind_width['overexpression'] == g[2]) &
        (ind_width['inducer_conc'] == g[3])]
    if (len(_ind_width) == 0):
        continue

    # find the x,y positions
    width10 = _ind_width[_ind_width['interval']
                         == '10%'][['lower', 'upper']].values.mean()
    frac10 = d[d['interval'] == '10%'][['lower', 'upper']].values.mean()
    ax.plot(frac10, width10, 'o', markeredgecolor='white',
            ms=4, markeredgewidth=0.5, color=cmaps[g[2]]['10%'], zorder=1000)
    for i, (_g, _d) in enumerate(d.groupby(['interval'], sort=False)):
        ax.hlines(width10, _d['lower'], _d['upper'], color=cmaps[g[2]][_g],
                  zorder=i + 1, label='__nolegend__', lw=0.1 + i * width_inc)

    for i, (_g, _d) in enumerate(_ind_width.groupby(['interval'], sort=False)):
        ax.vlines(frac10, _d['lower'], _d['upper'], color=cmaps[g[2]][_g],
                  zorder=i + 1, label='__nolegend__', lw=0.1 + i * width_inc)


ax.plot(mass_fracs['biomass_frac'], mass_fracs['width'], 'o')
phi_range = np.linspace(0.004, 0.1)
n_draws = 2000
uncertainty = np.zeros((n_draws, len(phi_range)))
for i in tqdm.tqdm(range(n_draws)):
    delta = np.random.normal(0.02, 0.001)
    alpha = 3
    # alpha = np.random.normal(3.2, 0.1)
    k = 0.1
    w_max = 1  # np.random.normal(0.9, 0.05)
    w_min = 0.25  # np.random.normal(0.25, 0.005)
    Lam = 12 * alpha * delta / (3 * alpha - 1)
    # phi_min = 0.01
    # phi_max = 0.08
    # phi_range = np.linspace(phi_min, phi_max)
    # phi_min = 0.01
    phi_min = Lam * k / (w_max + Lam * (k - 1))
    # phi_max = 0.08
    phi_max = Lam * k / (w_min + Lam * (k - 1))
    # del_phi = phi_range - phi_min
    # - Lam * k * del_phi**3 / phi_min**4 +  Lam * k * del_phi**5 / phi_min**6
    # pref = w_max - Lam * k * del_phi / phi_min**2

    # uncertainty[i, :] = pref
    slope = Lam * k / phi_max
    uncertainty[i, :] = w_max + slope * (1 - phi_range/phi_min)
    plt.plot(phi_range, uncertainty[i, :], 'k-', lw=0.1, alpha=0.1, zorder=1)
ax.plot([], [], '-', lw=1, color=cor['primary_green'], label='wildtype')
ax.plot([], [], '-', lw=1, color=cor['primary_blue'], label='malE-rbsB-fliC KO')
ax.plot([], [], '-', lw=1, color=cor['primary_gold'], label='rbsB OE in d3')
ax.plot([], [], '-', lw=1, color=cor['primary_purple'], label='malE OE in d3')
ax.plot([], [], '-', lw=1, color=cor['primary_black'], label='lacZ OE in WT')

ax.set_ylim([0.5, 1])
ax.set_xlim([0, 0.08])
ax.set_xlabel('$M_{peri} / M_{biomass}$')
ax.set_ylabel('width [µm]')

# ax.plot(phi_range, 1/pred,  'k-', lw=2)
ax.legend()
plt.savefig('/Users/gchure/Desktop/theory_fit_complete_analysis.pdf')


# %%
widths = samples.posterior.width_mu.to_dataframe().reset_index()
for g, d in size_data.groupby(['strain', 'carbon_source', 'overexpression', 'inducer_conc', 'size_cond_idx']):
    widths.loc[widths['width_mu_dim_0'] == g[-1] - 1, ['strain',
                                                       'carbon_source', 'overexpression', 'inducer_conc']] = g[:-1]


widths
