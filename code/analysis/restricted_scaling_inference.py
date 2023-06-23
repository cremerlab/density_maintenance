# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az
import size.viz
import size.fluxparity
import cmdstanpy
import scipy.stats
const = size.fluxparity.load_constants()
cor, pal = size.viz.matplotlib_style()

# Load the necessary data sets and restrict
si_data_ = pd.read_csv('../../data/literature/Si2017/si2017_size_phiRb.csv')
si_data = si_data_[si_data_['chlor_conc'] == 0].copy()
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
loc = 'inner membrane'
sa = {'inner membrane': 1, 'outer membrane': 1, 'membrane': 2}
mem = ms_data[ms_data['localization'] == loc]
mem.drop_duplicates(inplace=True)
peri = ms_data[ms_data['localization'] == 'periplasm']
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']
biomass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')

# Our data
data = pd.read_csv('../../data/mcmc/perturbation_parameter_percentiles.csv')
data_ = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
data = pd.concat([data, data_])
wt = data[(data['strain'] == 'wildtype') & (data['overexpression'] == 'none') &
          (data['inducer_conc'] == 0)]

# %%
width_popt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
size_data['mass_fraction'] = phiRb_popt[1] + \
    phiRb_popt[0] * size_data['growth_rate_hr']
phiRb_data['width_um'] = width_popt[1] + \
    width_popt[0] * phiRb_data['growth_rate_hr']
ms = ms_data[ms_data['localization'] == 'ribosomal sector'].copy()
ms.rename(columns={'dataset_name': 'source', 'width': 'width_um',
          'mass_frac': 'mass_fraction'}, inplace=True)
merged = pd.concat([size_data[['source', 'width_um', 'mass_fraction', 'growth_rate_hr']],
                   phiRb_data[['source', 'width_um',
                               'mass_fraction', 'growth_rate_hr']],
                   ms[['source', 'width_um', 'mass_fraction', 'growth_rate_hr']]])
merged = merged[(merged['source'] != 'Si et al. 2017') & (
    merged['source'] != 'Taheri-Araghi et al. 2015')]
size_data = size_data[size_data['source'] != 'Si et al. 2017']
size_data = size_data[size_data['source'] != 'Taheri-Araghi et al. 2015']
# %%
model = cmdstanpy.CmdStanModel(
    stan_file='./stan_models/restricted_scaling_inference.stan')

# %%
# Define PPc information
N_ppc = 500
nu_range = np.linspace(0, 25, N_ppc)
phiRb_range = size.fluxparity.phiRb_optimal_allocation(
    const['gamma_max'], nu_range, const['Kd_cpc'], const['phi_O'])
lam_range = size.fluxparity.steady_state_growth_rate(
    const['gamma_max'], phiRb_range, nu_range, const['Kd_cpc'], const['phi_O'])

# Define the data dict to condition the model on.
data_dict = {'N_size': len(merged),
             'N_alpha': len(size_data),
             'N_ms': len(mem),
             'N_biomass': len(biomass),
             'N_ppc': N_ppc,

             'width': merged['width_um'].values,
             'phiRb': merged['mass_fraction'].values,
             'alpha': size_data['aspect_ratio'].values,
             'phi_mem': mem['mass_frac'].values,
             'm_peri': peri['mass_fg'].values,
             'rho_biomass': biomass['drymass_density_fg_fL'].values,

             'beta_rp': 0.4558,
             'sa_prefactor': sa[loc],
             'delta': 0.0246,
             'phiRb_range': phiRb_range
             }
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(4, 3, figsize=(5, 6))
ax = ax.ravel()
# Use the bottom three rows for a big-ass legend
for a in ax[-3:]:
    a.axis('off')

# Set axis limits
ax[0].set_xlim([0, 0.7])
ax[0].set_ylim([0.4, 1.2])
ax[1].set_ylim([1, 6])
ax[2].set_ylim([0, 0.1])
ax[3].set_ylim([0.45, 1.2])
ax[4].set_ylim([1, 4])
ax[5].set_ylim([-0.1, 3])
ax[6].set_ylim([0, 4])
ax[7].set_ylim([0, 150])
ax[8].set_ylim([100, 400])

# Set axis labels
ax[0].set_xlabel('RNA / Protein', fontsize=6)
ax[0].set_ylabel('average cell width [µm]', fontsize=6)
ax[1].set_ylabel('average aspect ratio', fontsize=6)
ax[2].set_ylabel(f'allocation toward\n{loc} protein', fontsize=6)
ax[3].set_ylabel('average cell width [µm]', fontsize=6)
ax[4].set_ylabel('average cell length [µm]', fontsize=6)
ax[5].set_ylabel('average cell volume [µm$^3$]', fontsize=6)
ax[6].set_ylabel(f'{loc} protein density [fg / µm$^2$]', fontsize=6)
ax[7].set_ylabel(f'periplasmic protein density [fg / µm$^3$]', fontsize=6)
ax[8].set_ylabel(f'biomass density [fg / µm$^3$]', fontsize=6)
for a in ax[1:]:
    a.set_xlim([0, 2])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot the datasets and store the sources in a stupid way
sources = []

for g, d in merged.groupby('source'):
    if g not in sources:
        sources.append(g)
    _style = size.viz.style_point(g)
    ax[0].plot(d['mass_fraction'] / 0.4558, d['width_um'], **_style)

for g, d in size_data.groupby('source'):
    if g not in sources:
        sources.append(g)
    _style = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['aspect_ratio'], **_style)

for g, d in mem.groupby('dataset_name'):
    if g not in sources:
        sources.append(g)
    _style = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['mass_frac'], **_style)
    ax[6].plot(d['growth_rate_hr'], d['mass_fg'] /
               (sa[loc] * d['surface_area']), **_style)

for g, d in peri.groupby('dataset_name'):
    if g not in sources:
        sources.append(g)
    _style = size.viz.style_point(g)
    ax[7].plot(d['growth_rate_hr'], d['mass_fg'] /
               (sa[loc] * d['surface_area'] * 0.0246), **_style)

for g, d in size_data.groupby('source'):
    if g not in sources:
        sources.append(g)
    _style = size.viz.style_point(g)
    ax[3].plot(d['growth_rate_hr'], d['width_um'], **_style)
    ax[4].plot(d['growth_rate_hr'], d['length_um'], **_style)
    ax[5].plot(d['growth_rate_hr'], d['volume_um3'], **_style)

for g, d in biomass.groupby('source'):
    if g not in sources:
        sources.append(g)
    _style = size.viz.style_point(g)
    ax[8].plot(d['growth_rate_hr'], d['drymass_density_fg_fL'], **_style)

# Plot our data
pars = ['width_mu', 'length_mu', 'volume_mu',
        'alpha_mu', 'rho_peri', 'growth_mu']
axes = {'width_mu': ax[3], 'length_mu': ax[4],
        'volume_mu': ax[5], 'alpha_mu': ax[1],
        'rho_peri': ax[7]}
intervals = ['95%', '50%']
lam = wt[wt['quantity'] == 'growth_mu']
lam['phiRb_lower'] = (phiRb_popt[1] + phiRb_popt[0] * lam['lower']) / 0.4558
lam['phiRb_upper'] = (phiRb_popt[1] + phiRb_popt[0] * lam['upper']) / 0.4558
lam_med = lam[lam['interval'] == 'median']
for _p in pars[:-1]:
    if _p == 'rho_peri':
        lam = wt[wt['quantity'] == 'growth_mu']
        lam = lam[lam['carbon_source'] != 'LB']
        lam_med = lam[lam['interval'] == 'median']
    else:
        lam = wt[wt['quantity'] == 'growth_mu']
        lam_med = lam[lam['interval'] == 'median']

    med = wt[(wt['quantity'] == _p) & (wt['interval'] == 'median')]
    perc = wt[(wt['quantity'] == _p) & (wt['interval'] == '95%')]
    axes[_p].vlines(lam_med['lower'], perc['lower'],
                    perc['upper'], lw=1, color=cor['primary_blue'])
    axes[_p].hlines(med['lower'], lam[lam['interval'] == '95%']['lower'], lam[lam['interval'] == '95%']['upper'],
                    lw=1, color=cor['primary_blue'])

    axes[_p].plot(lam_med['lower'], med['lower'], 'o', ms=3, markeredgecolor=cor['primary_blue'],
                  color='white', markeredgewidth=1)

lam = wt[wt['quantity'] == 'growth_mu']
lam['phiRb_lower'] = (phiRb_popt[1] + phiRb_popt[0] * lam['lower']) / 0.4558
lam['phiRb_upper'] = (phiRb_popt[1] + phiRb_popt[0] * lam['upper']) / 0.4558
lam_med = lam[lam['interval'] == 'median']

med = wt[(wt['quantity'] == 'width_mu') & (wt['interval'] == 'median')]
ax[0].vlines(lam_med['phiRb_lower'], med['lower'],
             med['upper'], lw=1, color=cor['primary_blue'])
ax[0].hlines(med['lower'], lam[lam['interval'] == '95%']['phiRb_lower'],
             lam[lam['interval'] == '95%']['phiRb_upper'], lw=1, color=cor['primary_blue'])
ax[0].plot(lam_med['phiRb_lower'], med['lower'], 'o', ms=3, markeredgecolor=cor['primary_blue'],
           color='white', markeredgewidth=1)

# Plot the fits and predictions
percs = [(2.5, 97.5), (12.5, 87.5), (37.5,  62.5), (50.0, 50.0)]
hue = ['pale', 'light', 'primary', 'dark']
color = ['black', 'blue', 'blue', 'blue']
ppcs = ['width_rep', 'ell_rep', 'vol_rep', 'rho_peri_rep']
axes = {'width_rep': [ax[0], ax[3]], 'ell_rep': [
    ax[4]], 'vol_rep': [ax[5]], 'rho_peri_rep': [ax[7]]}
df = samples.posterior[['alpha_rep', 'phi_mem_rep', 'rho_mem', 'rho_biomass_rep']
                       ].to_dataframe().reset_index()
df.dropna(inplace=True)
for i, v in enumerate(['alpha_rep', 'phi_mem_rep', 'rho_mem', 'rho_biomass_rep']):
    if i == 2:
        _ = 4
        c = 'blue'
    elif i == 3:
        _ = 5
        c = 'black'
    else:
        _ = 1
        c = 'black'
    for j, p in enumerate(percs):
        lower, upper = np.percentile(df[v], p)
        ax[i + _].fill_between(lam_range, lower, upper,
                               color=cor[f'{hue[j]}_{c}'], alpha=0.5)

for i, par in enumerate(ppcs):
    df = samples.posterior[par].to_dataframe().reset_index()
    df.dropna(inplace=True)
    for j, p in enumerate(percs):
        lower, upper = [], []
        for g, d in df.groupby(f'{par}_dim_0'):
            _p = np.percentile(d[par], p)
            lower.append(_p[0])
            upper.append(_p[1])
        if i == 0:
            x = phiRb_range / 0.4558
            _ax = axes[par][0]
            axes[par][1].fill_between(lam_range, lower, upper, color=cor[f'{hue[j]}_blue'],
                                      alpha=0.5)
        else:
            x = lam_range
            _ax = axes[par][0]
        _ax.fill_between(x, lower, upper, color=cor[f'{hue[j]}_{color[i]}'],
                         alpha=0.5)


for s in sources:
    _style = size.viz.style_point(s)
    ax[9].plot([], [], ms=4, **_style)
ax[9].plot([], [], 'o', ms=3,
           markeredgecolor=cor['primary_blue'], color='w', markeredgewidth=1, label='this study')
ax[9].legend(fontsize=4, ncols=5, handletextpad=0.1, bbox_to_anchor=(3.5, 0.5))

plt.tight_layout()
plt.savefig(
    f'../../figures/FigX_{loc}_restricted_inference.pdf', bbox_inches='tight')
# %%
shift_data = pd.read_csv('../../data/literature/woldringh1980_upshift.csv')
# %%
# Shift dynamics
lam_pre = np.log(2) / (72/60)  # Preshift with doubling time of 72 min
lam_post = np.log(2) / (24/60)  # Postshift with doubling time of 24 min

# Estimate the phiRb given linear relation
phiRb_pre = phiRb_popt[1] + phiRb_popt[0] * lam_pre
phiRb_post = phiRb_popt[1] + phiRb_popt[0] * lam_post

# Estimate nu max based on this
nu_pre = size.fluxparity.estimate_nu_FPM(
    phiRb_pre, lam_pre, const, const['phi_O'])
nu_post = size.fluxparity.estimate_nu_FPM(
    phiRb_post, lam_post, const, const['phi_O'])

# %%
preshift = {'gamma_max': const['gamma_max'],
            'nu_max': nu_pre,
            'Kd_TAA': const['Kd_TAA'],
            'Kd_TAA_star': const['Kd_TAA_star'],
            'kappa_max': const['kappa_max'],
            'phi_O': const['phi_O'],
            'tau': const['tau']}
postshift = {'gamma_max': const['gamma_max'],
             'nu_max': nu_post,
             'Kd_TAA': const['Kd_TAA'],
             'Kd_TAA_star': const['Kd_TAA_star'],
             'kappa_max': const['kappa_max'],
             'phi_O': const['phi_O'],
             'tau': const['tau']}

df = size.fluxparity.nutrient_shift_FPM([preshift, postshift])

df['Rb_content'] = df['M_Rb'] / df['M']
alpha = float(np.mean(samples.posterior.alpha_mu).values)
k = float(np.mean(samples.posterior.k).values)
phi_mem = float(np.mean(samples.posterior.phi_mem_mu).values)
pref = (12 * alpha * sa[loc]) / (k * phi_mem * (3 * alpha - 1))
df['width'] = pref * (1 + df['Rb_content'].values / 0.4558)
# %%
fig, ax = plt.subplots(2, 1, figsize=(4, 3), sharex=True)
ax[0].set_ylim([0, 0.25])
ax[1].set_ylim([0.5, 1])
ax[0].plot(df['shifted_time'], df['Rb_content'],
           '--', lw=1, color=cor['primary_red'], label='flux-parity parity prediction')
ax[1].plot(df['shifted_time'], df['width'],
           '--', lw=1, color=cor['primary_red'], label='density maintenance prediction')
ax[1].plot(shift_data['time_hr'], shift_data['width_um'], 'o', ms=4,
           markeredgecolor=cor['primary_black'], color=cor['light_black'], alpha=0.5,
           markeredgewidth=0.25, label='Woldringh 1980')

ax[0].legend(fontsize=6)
ax[1].legend(fontsize=6)
ax[1].set_xlabel('time from upshift [hr]', fontsize=6)
ax[0].set_ylabel('ribosome content', fontsize=6)
ax[1].set_ylabel('cell width [µm]', fontsize=6)

# %%
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
ax.set_xlim([0, 2])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('RNA / Protein', fontsize=6)

for g, d in merged.groupby('source'):
    _style = size.viz.style_point(g, alpha=0.4)
    ax.plot(d['growth_rate_hr'], d['mass_fraction']/0.4558, ms=4.5, **_style)

ax.plot(lam_range, phiRb_range/0.4558, '-', lw=1,
        color=cor['primary_red'], label='optimal allocation')
ax.plot(lam_range, (phiRb_popt[1] + phiRb_popt[0] * lam_range) / 0.4558,
        '-', lw=1, color=cor['primary_blue'], label='empirical relation')

ax.legend(fontsize=4, bbox_to_anchor=(1, 1))
plt.savefig('../../figures/FigXG_RNA_P_scaling.pdf')

# %%
