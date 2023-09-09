# %%
import numpy as np
import scipy.stats
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import size.fluxparity
cor, pal = size.viz.matplotlib_style()

peri_width = pd.read_csv(
    '../../data/literature/Asmar2017/asmar2017_wt_periplasm_width.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Si et al. 2017']
prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
vol_data = pd.read_csv(
    '../../data/literature/collated_literature_volume_data.csv')
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
ms_emp = pd.read_csv(
    '../../data/mcmc/mass_spec_empirical_summaries_wide.csv')

# Compute empirical fits
phiRb_popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
vol_popt = scipy.stats.linregress(
    vol_data['growth_rate_hr'], np.log(vol_data['volume_um3']))

lam_range = np.linspace(0, 2.5, 200)
phiRb_fit = phiRb_popt[1] + phiRb_popt[0] * lam_range
vol_fit = np.exp(vol_popt[1] + vol_popt[0] * lam_range)

# Set up the figure canvas for the "growth laws"
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
ax[-1].axis('off')
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]\n$\lambda$', fontsize=6)
ax[0].set_ylabel('$\phi_{Rb}$\nribosomal protein mass fraction', fontsize=6)
ax[1].set_ylabel('$V$\naverage cell volume [µm$^3$]', fontsize=6)


# force limits
ax[0].set_ylim([0, 0.3])

# Plot the phiRb data
labels = []
for g, d in phiRb_data.groupby('source'):
    if g not in labels:
        labels.append(g)
    fmt = size.viz.style_point(g, alpha=1)
    ax[0].plot(d['growth_rate_hr'], d['mass_fraction'], **fmt)

ax[0].plot(lam_range, phiRb_fit, '--', color=cor['light_blue'], lw=2)

# Plot the size data
for g, d in vol_data.groupby('source'):
    if g not in labels:
        labels.append(g)
    fmt = size.viz.style_point(g, alpha=1)
    ax[1].plot(d['growth_rate_hr'], np.log2(d['volume_um3']), **fmt)

ax[1].plot(lam_range, np.log2(vol_fit), '--', color=cor['light_blue'], lw=2)

for l in labels:
    fmt = size.viz.style_point(l)
    ax[2].plot([], [], ms=4, **fmt)

ax[2].legend(fontsize=5, handletextpad=0.5)
ticks = [-2, -1, 0, 1, 2]
ax[1].set_ylim([-2, 2.5])
ax[1].set_yticks(ticks)
ax[1].set_yticklabels([f'{2**t:0.1f}' for t in ticks])
plt.tight_layout()
plt.savefig('../../figures/Fig1_growth_law_plots.pdf', bbox_inches='tight')


# %%
# Plot the empirical constants from
# Note that this trims the displayed data from a growth rate of 0 - 1.5. A few
# pints lay past this, but it squishes the data on the low end, making trends
# difficult to see. Will include the full range as a SI fig


fig, ax = plt.subplots(2, 2, figsize=(4, 3))
ax = ax.ravel()
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]\n$\lambda$', fontsize=6)
    a.set_xlim([0, 1.5])
ax[0].set_ylabel('$\phi_{mem}$\nmembrane protein\nmass fraction', fontsize=6)
ax[1].set_ylabel(
    r'$\rho_{mem}$'+'\nmembrane protein\nareal density [fg / µm$^2$]', fontsize=6)
ax[2].set_ylabel(
    '$\phi_{peri}$\nperiplasmic protein\nmass fraction', fontsize=6)
ax[3].set_ylabel(
    r'$\rho_{peri}$'+'\nperiplasmic\nprotein density [fg / µm$^3$]', fontsize=6)

ax[0].set_ylim([0, 0.25])
ax[1].set_ylim([0, 5])
ax[2].set_ylim([0, 0.1])
# ax[3].set_ylim([0, 20])


for g, d in ms_data.groupby('dataset_name'):
    mem = d[d['localization'] == 'membrane']
    peri = d[d['localization'] == 'periplasm']
    fmt = size.viz.style_point(g)
    fmt['color'] = cor['light_blue']
    ax[0].plot(mem['growth_rate_hr'], mem['mass_frac'], **fmt)
    fmt['color'] = cor['light_purple']
    ax[2].plot(peri['growth_rate_hr'], peri['mass_frac'], **fmt)

for g, d in ms_emp.groupby('source'):
    rho_mem = d[d['quantity'] == 'ms_rho_mem']
    rho_peri = d[d['quantity'] == 'ms_rho_peri']
    fmt = size.viz.style_point(g)
    fmt['label'] = '__nolegend__'

    fmt['color'] = cor['light_blue']
    ax[1].vlines(rho_mem['growth_rate_hr'], rho_mem['2.5%'],
                 rho_mem['97.5%'], color=fmt['color'], linewidth=0.5)
    ax[1].plot(rho_mem['growth_rate_hr'], rho_mem['median_value'], **fmt)

    fmt['color'] = cor['light_purple']
    ax[3].vlines(rho_peri['growth_rate_hr'], rho_peri['2.5%'],
                 rho_peri['97.5%'], color=fmt['color'], linewidth=0.5)
    ax[3].plot(rho_peri['growth_rate_hr'], rho_peri['median_value'], **fmt)
    fmt['label'] = g
    ax[1].plot([], [], ms=4, **fmt)
ax[1].legend(fontsize=5, handletextpad=0.5)
# Plot the trends

# Fit the periplasmic linear regime
peri = ms_data[(ms_data['localization'] == 'periplasm')
               & (ms_data['growth_rate_hr'] <= 1.5)]
peri_popt = scipy.stats.linregress(peri['growth_rate_hr'], peri['mass_frac'])
peri_emp = ms_emp[(ms_emp['quantity'] == 'ms_rho_peri')
                  & (ms_emp['growth_rate_hr'] <= 1)]
rho_peri_popt = scipy.stats.linregress(
    peri_emp['growth_rate_hr'], peri_emp['median_value'])
lam_range = np.linspace(0, 1.5)
peri_fit = peri_popt[1] + peri_popt[0] * lam_range
rho_peri_fit = rho_peri_popt[1] + rho_peri_popt[0] * lam_range

ax[0].hlines(0.12, 0, 1.5, color=cor['black'],
             linestyle='--', linewidth=2)
ax[1].hlines(2.5, 0, 1.5, color=cor['black'], linestyle='--', linewidth=2)
ax[2].plot(lam_range, peri_fit, color=cor['black'],
           linestyle='--', linewidth=2)
ax[3].plot(lam_range, rho_peri_fit, color=cor['black'],
           linewidth=2, linestyle='--')
# ax[3].hlines(10, 0, 1.5, color=cor['light_blue'], linestyle='--', linewidth=2)

# plt.tight_layout()
plt.savefig('../../figures/Fig1_empirical_constants.pdf', bbox_inches='tight')

# %%

prot_popt = scipy.stats.linregress(
    prot_data['growth_rate_hr'], np.log(prot_data['fg_protein_per_cell']))
_size_data = size_data[size_data['growth_rate_hr'] <= 2]
sa_popt = scipy.stats.linregress(
    _size_data['growth_rate_hr'], _size_data['surface_area_um2'])
lam_range = np.linspace(0, 2.5, 200)
prot_fit = np.exp(prot_popt[1] + prot_popt[0] * lam_range)
sa_fit = sa_popt[1] + sa_popt[0] * lam_range

fig, ax = plt.subplots(2, 1, figsize=(1, 2), sharex=True)
ax[1].set_ylim([50, 800])
ax[0].set_xlim([0, 2])
ax[1].set_xlim([0, 2])
ax[0].set_ylim([0, 13])
ax[0].set_ylabel('$S_A$\nsurface area [µm$^2$]', fontsize=6)
ax[1].set_ylabel('$M_{prot}$\ntotal protein [fg/cell]', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

# ax[2].set_yticks([])
# ax[2].set_ylabel('frequency', fontsize=6)
# ax[2].set_xlabel('periplasm width [nm]\n$\delta$', fontsize=6)

for g, d in prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[1].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt, ms=4)
ax[1].plot(lam_range, prot_fit, '-', lw=1, color=cor['primary_blue'])
for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt, ms=4)
ax[0].plot(lam_range, sa_fit, '-', lw=1, color=cor['primary_blue'])

# ax[2].bar(peri_width['width_nm'], peri_width['frequency'],
#   color=cor['light_black'])
# ax[2].vlines(24.6, 0, 0.2, color=cor['primary_blue'], linewidth=1)
ax[0].legend(fontsize=6)
plt.savefig('../../figures/Fig1_empirical_relations.pdf', bbox_inches='tight')
