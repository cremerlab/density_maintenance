# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az
import size.viz
import size.fluxparity
import cmdstanpy
const = size.fluxparity.load_constants()
cor, pal = size.viz.matplotlib_style()

# Load the necessary data sets and restrict
si_data_ = pd.read_csv('../../data/literature/Si2017/si2017_size_phiRb.csv')
si_data = si_data_[si_data_['chlor_conc'] == 0].copy()
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
mem = ms_data[ms_data['localization'] == 'membrane']
mem.drop_duplicates(inplace=True)
peri = ms_data[ms_data['localization'] == 'periplasm']
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data['aspect_ratio'] = size_data['length_um'] / size_data['width_um']
biomass = pd.read_csv('../../data/literature/collated_drymass_densities.csv')


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
data_dict = {'N_size': len(si_data),
             'N_alpha': len(size_data),
             'N_ms': len(mem),
             'N_biomass': len(biomass),
             'N_ppc': N_ppc,

             'width': si_data['width_um'].values,
             'phiRb': si_data['phi_Rb'].values,
             'alpha': size_data['aspect_ratio'].values,
             'phi_mem': mem['mass_frac'].values,
             'm_peri': peri['mass_fg'].values,
             'rho_biomass': biomass['drymass_density_fg_fL'].values,

             'beta_rp': 0.4558,
             'sa_prefactor': 2,
             'delta': 0.0246,
             'phiRb_range': phiRb_range
             }
_samples = model.sample(data=data_dict, show_console=True)
samples = az.from_cmdstanpy(_samples)

# %%
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
ax = ax.ravel()

# Set axis limits
ax[0].set_xlim([0, 0.3])
ax[0].set_ylim([0.4, 1.2])
ax[1].set_ylim([1, 10])
ax[2].set_ylim([0, 0.15])
ax[3].set_ylim([0.45, 1.2])
ax[4].set_ylim([1, 6])
ax[5].set_ylim([-0.1, 5])

# Set axis labels
ax[0].set_xlabel('allocation towards ribosomes', fontsize=6)
ax[0].set_ylabel('average cell width [µm]', fontsize=6)
ax[1].set_ylabel('average aspect ratio', fontsize=6)
ax[2].set_ylabel('allocation toward\ninner membrane protein', fontsize=6)
ax[3].set_ylabel('average cell width [µm]', fontsize=6)
ax[4].set_ylabel('average cell length [µm]', fontsize=6)
ax[5].set_ylabel('average cell volume [µm$^3$]', fontsize=6)
for a in ax[1:]:
    a.set_xlim([0, 2.5])
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)


# Plot data
ax[0].plot(si_data['phi_Rb'], si_data['width_um'],
           **size.viz.style_point('Si et al. 2017'))

for g, d in size_data.groupby('source'):
    _style = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['aspect_ratio'], **_style)

for g, d in mem.groupby('dataset_name'):
    _style = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['mass_frac'], **_style)

for g, d in size_data.groupby('source'):
    _style = size.viz.style_point(g)
    ax[3].plot(d['growth_rate_hr'], d['width_um'], **_style)
    ax[4].plot(d['growth_rate_hr'], d['length_um'], **_style)
    ax[5].plot(d['growth_rate_hr'], d['volume_um3'], **_style)


# Plot the fits and predictions
percs = [(2.5, 97.5), (12.5, 87.5), (37.5,  62.5), (50.0, 50.0)]
hue = ['pale', 'light', 'primary', 'dark']
color = ['black', 'blue', 'blue', 'blue']
ppcs = ['width_rep', 'ell_rep', 'vol_rep']
axes = {'width_rep': [ax[0], ax[3]], 'ell_rep': [ax[4]], 'vol_rep': [ax[5]]}

df = samples.posterior[['alpha_rep', 'phi_mem_rep']
                       ].to_dataframe().reset_index()
df.dropna(inplace=True)
for i, v in enumerate(['alpha_rep', 'phi_mem_rep']):
    for j, p in enumerate(percs):
        lower, upper = np.percentile(df[v], p)
        ax[i + 1].fill_between(lam_range, lower, upper,
                               color=cor[f'{hue[j]}_black'], alpha=0.5)

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
            x = phiRb_range
            ax = axes[par][0]
            axes[par][1].fill_between(lam_range, lower, upper, color=cor[f'{hue[j]}_blue'],
                                      alpha=0.5)
        else:
            x = lam_range
            ax = axes[par][0]
        ax.fill_between(x, lower, upper, color=cor[f'{hue[j]}_{color[i]}'],
                        alpha=0.5)

plt.tight_layout()
