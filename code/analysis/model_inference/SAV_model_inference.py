# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./SAV_model_inference.stan')

size_data = pd.read_csv(
    '../../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Si et al. 2017']
prot_data = pd.read_csv(
    '../../../data/literature/collated_protein_per_cell.csv')
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')
mem_data = ms_data[ms_data['localization'] == 'membrane']
peri_data = ms_data[ms_data['localization'] == 'periplasm']

phiRb_data = pd.read_csv(
    '../../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']

density_data = pd.read_csv(
    '../../../data/literature/collated_drymass_densities.csv')


N_pred = 300
pred_phiRb = np.linspace(0.05, 0.30, N_pred)
data_dict = {
    'N_pred': N_pred,
    'pred_phiRb': pred_phiRb,

    'N_size': len(size_data),
    'size_lam': size_data['growth_rate_hr'].values,
    'surface_areas': size_data['surface_area_um2'].values,

    'N_prot': len(prot_data),
    'prot_lam': prot_data['growth_rate_hr'].values,
    'prot_per_cell': prot_data['fg_protein_per_cell'].values,

    'N_drymass': len(density_data),
    'drymass_lam': density_data['growth_rate_hr'].values,
    'drymass_density': density_data['drymass_density_fg_fL'].values,

    'N_ms': len(mem_data),
    'ms_lam': mem_data['growth_rate_hr'].values,
    'phi_mem': mem_data['mass_frac'].values,
    'phi_peri':  peri_data['mass_frac'].values,

    'N_phiRb': len(phiRb_data),
    'phiRb': phiRb_data['mass_fraction'].values,
    'phiRb_lam': phiRb_data['growth_rate_hr'].values,
}

_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
wt = pd.read_csv(
    '../../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
# %%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_ylabel('surface to volume [µm$^{-1}$]', fontsize=6)
ax.set_xlabel('ribosomal allocation', fontsize=6)
labels = []
for g, d in samples.posterior.size_phiRb.to_dataframe().reset_index().groupby('size_phiRb_dim_0'):
    src = size_data['source'].values[g]
    if src not in labels:
        labels.append(src)
        lab = src
    else:
        lab = '__nolegend__'
    fmt = size.viz.style_point(size_data['source'].values[g])
    fmt['label'] = lab
    med = np.median(d['size_phiRb'].values)
    perc = np.percentile(d['size_phiRb'].values, (2.5, 97.5))
    ax.hlines(size_data['surface_to_volume'].values[g], perc[0],
              perc[1], linewidth=1, color=cor['primary_black'])
    ax.plot(med, size_data['surface_to_volume'].values[g], ms=4, **fmt)
ax.plot([], [], 'o', markeredgecolor=cor['dark_blue'], markerfacecolor='w',
        ms=4, markeredgewidth=1, label='this study')
ax.legend(fontsize=4)
lowers = np.array([2.5, 12.5, 62.5, 50])
uppers = 100 - lowers
labels = ['95%', '75%', '25%', 'median']
theory_post = samples.posterior.theory.to_dataframe().reset_index()
percs = size.viz.compute_percentiles(theory_post, 'theory', 'theory_dim_0',
                                     lower_bounds=lowers, upper_bounds=uppers,
                                     interval_labels=labels)

int_colors = {'95%': cor['pale_blue'],
              '75%': cor['light_blue'],
              '25%': cor['primary_blue'],
              'median': cor['blue']}
for g, d in percs.groupby('interval', sort=False):
    if g != 'median':
        ax.fill_between(pred_phiRb, d['lower'],
                        d['upper'], color=int_colors[g], alpha=0.5)
    else:
        ax.plot(pred_phiRb, d['lower'],
                color=int_colors[g], linewidth=1, alpha=0.5)

for g, d in wt.groupby('carbon_source'):
    phiRb = d[d['quantity'] == 'phi_Rb']
    sav = d[d['quantity'] == 'surface_to_volume']
    ax.vlines(phiRb['median_value'], sav['2.5%'], sav['97.5%'],
              linewidth=1, color=cor['blue'])
    ax.hlines(sav['median_value'], phiRb['2.5%'], phiRb['97.5%'],
              linewidth=1, color=cor['blue'])
    ax.plot(phiRb['median_value'], sav['median_value'], 'o', markeredgecolor=cor['blue'],
            markerfacecolor='w', markeredgewidth=1, ms=4)
