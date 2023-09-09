# %%
import size.fluxparity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import scipy.stats
cor, pal = size.viz.matplotlib_style()

si_data = pd.read_csv('../../data/literature/Si2017/si2017_chlor_phiRb.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
wpopt = scipy.stats.linregress(
    size_data['growth_rate_hr'], size_data['width_um'])
phiRb_data['width'] = wpopt[1] + wpopt[0] * phiRb_data['growth_rate_hr']
wt_params = pd.read_csv(
    '../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
pert_params = pd.read_csv(
    '../../data/mcmc/perturbation_posterior_parameter_summaries.csv')
buke = pd.read_csv('../../data/literature/Buke2022/buke2022_processed.csv')

fig, ax = plt.subplots(3, 2, figsize=(6, 4))
ax = ax.ravel()
ax_mapper = {'phi_mem': ax[0], 'phi_peri': ax[1],
             'rho_mem': ax[2], 'rho_peri': ax[3],
             'phi_Rb': ax[4]}
ax[0].set_ylim([0, 0.25])
ax[1].set_ylim([0, 0.15])
ax[2].set_ylim([0, 5])
ax[3].set_ylim([0, 250])
ax[4].set_ylim([0, 0.35])

for a in ax[:-1]:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('membrane allocation', fontsize=6)
ax[1].set_ylabel('periplasm allocation', fontsize=6)
ax[2].set_ylabel('membrane density [fg / µm$^2$]', fontsize=6)
ax[3].set_ylabel('periplasm density [fg / µm$^3$]', fontsize=6)
ax[4].set_ylabel('ribosomal allocation', fontsize=6)
ax[5].set_xlabel('ribosomal allocation', fontsize=6)
ax[5].set_ylabel('width [µm]', fontsize=6)

# Plot the literature data
for g, d in ms_data.groupby('dataset_name'):
    mem = d[d['localization'] == 'membrane']
    mem['rho'] = mem['mass_fg'] / (2 * mem['surface_area'])
    peri = d[d['localization'] == 'periplasm']
    peri['rho'] = peri['mass_fg'] / (peri['surface_area'] * 0.0246)
    fmt = size.viz.style_point(g, alpha=0.25)
    ax[0].plot(mem['growth_rate_hr'], mem['mass_frac'], **fmt)
    ax[1].plot(peri['growth_rate_hr'], peri['mass_frac'], **fmt)
    ax[2].plot(mem['growth_rate_hr'], mem['rho'], **fmt)
    ax[3].plot(peri['growth_rate_hr'], peri['rho'], **fmt)

for g, d in phiRb_data.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.25)
    ax[4].plot(d['growth_rate_hr'], d['mass_fraction'], **fmt)
    ax[5].plot(d['mass_fraction'], d['width'], **fmt)

wt_lam = wt_params[wt_params['quantity'] == 'growth_rate']
for i, (g, d) in enumerate(wt_params[
        wt_params['quantity'].isin(['phi_mem', 'phi_peri',
                                    'rho_mem', 'rho_peri',
                                    'phi_Rb'])].groupby(['quantity'])):

    lam = wt_lam.copy()
    if 'peri' in g:
        lam = lam[lam['carbon_source'] != 'LB']
    ax_mapper[g].hlines(d['mean_value'], lam['2.5%'], lam['97.5%'], linewidth=1,
                        color=cor['primary_blue'])
    ax_mapper[g].vlines(lam['mean_value'], d['2.5%'], d['97.5%'], linewidth=1,
                        color=cor['primary_blue'])
    ax_mapper[g].plot(lam['mean_value'], d['mean_value'], 'o', markeredgewidth=1, markerfacecolor='w',
                      markeredgecolor=cor['primary_blue'])


wt_phiRb = wt_params[wt_params['quantity'] == 'phi_Rb']
wt_width = wt_params[wt_params['quantity'] == 'width']
ax[-1].vlines(wt_phiRb['mean_value'], wt_width['2.5%'], wt_width['97.5%'], linewidth=1,
              color=cor['primary_blue'])
ax[-1].hlines(wt_width['mean_value'], wt_phiRb['2.5%'], wt_phiRb['97.5%'], linewidth=1,
              color=cor['primary_blue'])
ax[-1].plot(wt_phiRb['mean_value'], wt_width['mean_value'], 'o', markeredgewidth=1,
            markeredgecolor=cor['primary_blue'], markerfacecolor='w')

pert_lam = pert_params[pert_params['quantity'] == 'growth_rate']
for i, (g, d) in enumerate(pert_params[
        pert_params['quantity'].isin(['phi_mem', 'phi_peri',
                                      'rho_mem', 'rho_peri',
                                      'phi_Rb'])].groupby(['quantity', 'overexpression'])):
    lam = pert_lam[pert_lam['overexpression'] == g[1]].copy()
    if g[1] == 'relA':
        c = cor['primary_red']
    else:
        c = cor['primary_gold']
    ax_mapper[g[0]].hlines(d['mean_value'], lam['2.5%'], lam['97.5%'], linewidth=1,
                           color=c)
    ax_mapper[g[0]].vlines(lam['mean_value'], d['2.5%'], d['97.5%'], linewidth=1,
                           color=c)
    ax_mapper[g[0]].plot(lam['mean_value'], d['mean_value'], 'o', markeredgewidth=1, markerfacecolor='w',
                         markeredgecolor=c)

for p, c in zip(['relA', 'meshI'], [cor['primary_red'], cor['primary_gold']]):
    p_phiRb = pert_params[(pert_params['quantity'] == 'phi_Rb') &
                          (pert_params['overexpression'] == p)]
    p_width = pert_params[(pert_params['quantity'] == 'width') &
                          (pert_params['overexpression'] == p)]
    ax[-1].vlines(p_phiRb['mean_value'], p_width['2.5%'], p_width['97.5%'], linewidth=1,
                  color=c)
    ax[-1].hlines(p_width['mean_value'], p_phiRb['2.5%'], p_phiRb['97.5%'], linewidth=1,
                  color=c)
    ax[-1].plot(p_phiRb['mean_value'], p_width['mean_value'], 'o', markeredgewidth=1,
                markeredgecolor=c, markerfacecolor='w')


# _si_data = si_data[si_data['carbon_source']
#                    == 'MOPS glucose'].groupby('chlor_conc').mean()
# __si_data = si_data[si_data['chlor_conc'] == 0]
# ax[-1].plot(_si_data['phi_Rb'], _si_data['width'],
# 'X', color=cor['primary_black'])
#
# ax[-2].plot(__si_data['growth_rate_hr'], __si_data['phi_Rb'],
# 'X', color=cor['primary_black'])

# Buke wt
# lam_wt = 1
# phiRb_wt = 0.12
# w_wt = 0.8

# lam_100 = 0.5
# phiRb_100 = 0.14
# w_100 = 1.1 * w_wt
# ax[-2].plot(lam_wt, phiRb_wt, 'X', color=cor['primary_green'])
# ax[-1].plot(phiRb_wt, w_wt, 'X', color=cor['primary_green'])
# ax[-2].plot(lam_100, phiRb_100, 'X', color=cor['primary_red'])
# ax[-1].plot(phiRb_100, w_100, 'X', color=cor['primary_red'])
# buke = buke[(buke['overexpression'] == 'meshI') & (buke['inducer_conc'] == 100)]
buke['color'] = [cor['primary_red'], cor['primary_black']]
for g, d in buke.groupby(['inducer_conc']):
    ax[-1].plot(d['phiRb'], d['width'], 'X', color=d['color'].values[0])
    ax[-2].plot(np.log(2) * d['lam'], d['phiRb'],
                'X', color=d['color'].values[0])

plt.tight_layout()


# plt.savefig('./perturbation_assessment.pdf', bbox_inches='tight')
