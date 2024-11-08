#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz
import scipy.stats
cor, pal = size.viz.matplotlib_style()

# Load the literature size data
lit_size_data = pd.read_csv('../../data/literature/full_literature_size_data.csv')
lit_size_data = lit_size_data[lit_size_data['source'] != 'Basan et al. 2015']
lit_prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')


# Load our coalesced data
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']=='wildtype']
prot_data = pd.read_csv('../../data/bulk_protein_per_cell.csv')

fig, ax = plt.subplots(1, 3, figsize=(6, 2))

for g, d in lit_size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['volume_um3'], **fmt)
    ax[1].plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt)

for g, d in lit_prot_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[2].plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

ax[2].errorbar(prot_data['mean_growth_rate_hr'], prot_data['fg_prot_per_cell'], 
                yerr=prot_data['err_fg_prot_per_cell'], 
                xerr=prot_data['std_growth_rate_hr'], 
                marker='D', markerfacecolor='w', markeredgecolor=cor['primary_black'],
                markeredgewidth=0.75, markersize=4, linewidth=0.5, color=cor['primary_black'],
                linestyle='none')

fmt = size.viz.style_point('This Study')
ax[0].errorbar(data['growth_rate_hr'], data['volume'], 
                xerr=data['growth_rate_hr_std'],
                    **fmt)
ax[1].errorbar(data['growth_rate_hr'], data['surface_area'], 
                xerr=data['growth_rate_hr_std'],
                **fmt)

# Determine empirical fits
lam_range = np.linspace(0, 2.5, 200)
vol_popt = scipy.stats.linregress(data['growth_rate_hr'], np.log(data['volume']))
vol_fit = np.exp(vol_popt.intercept + vol_popt.slope * lam_range)
sa_popt = scipy.stats.linregress(data['growth_rate_hr'], np.log(data['surface_area']))
sa_fit = np.exp(sa_popt.intercept + sa_popt.slope * lam_range)
prot_fit = scipy.stats.linregress(prot_data['mean_growth_rate_hr'], np.log(prot_data['fg_prot_per_cell']))
prot_fit = np.exp(prot_fit.intercept + prot_fit.slope * lam_range)

# Plot the fits
ax[0].plot(lam_range, vol_fit, '--', lw=1, color=cor['primary_black'], zorder=10)
ax[1].plot(lam_range, sa_fit, '--', lw=1, color=cor['primary_black'], zorder=10)
ax[2].plot(lam_range, prot_fit, '--', lw=1, color=cor['primary_black'], zorder=10)
ax[0].set_ylim([-0.5, 5])
# ax[1].set_ylim([1, 10])
ax[2].set_ylim([50, 1000])

# Add labels
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('volume [µm$^{3}$]', fontsize=6)
ax[1].set_ylabel('surface area [µm$^{2}$]', fontsize=6)
ax[2].set_ylabel('protein [fg/cell]', fontsize=6)
ax[0].legend()
ax[2].legend()
plt.savefig('./plots/fig2_size_relations.pdf', bbox_inches='tight')


#%%
# Load the estimated densities
densities = pd.read_csv('../analysis/output/mass_spec_densities_summary.csv')
densities = densities[densities['interval'].isin(['median', '95%'])]
fig, ax = plt.subplots(1, 3, figsize=(6, 1.5))
quantities = ['rho_cyt', 'sigma_mem', 'rho_peri']
colors = ['black', 'blue', 'purple']
for g, d in densities.groupby('source'):
    for i, q in enumerate(quantities):
        _d = d[d['quantity']==q]
        med = _d[_d['interval']=='median']
        perc = _d[_d['interval']=='95%']
        fmt = size.viz.style_point(g, alpha=0.4)
        if g == 'This Study':
            fmt['color'] = cor[f'pale_{colors[i]}']
            fmt['markeredgecolor'] = cor[f'primary_{colors[i]}']
            alpha = 1
        else:
            fmt['color'] = cor[f'primary_{colors[i]}']
            fmt['markeredgecolor'] = cor[colors[i]]
            alpha = 0.5
        ax[i].vlines(med['growth_rate_hr'], perc['lower'], perc['upper'], 
                lw=fmt['markeredgewidth'], alpha=alpha, color=fmt['markeredgecolor'])
        ax[i].plot(med['growth_rate_hr'], med['lower'], **fmt)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylim([100, 700])
ax[1].set_ylim([0, 5])
ax[2].set_ylim([0, 250])
ax[0].set_ylabel(r'$\rho_{cyt}$' + '\ncytoplasmic density [fg/fL]', fontsize=6)
ax[1].set_ylabel(r'$\sigma_{mem}$' + '\nmembrane density [fg/µm$^2$]', fontsize=6)
ax[2].set_ylabel(r'$\rho_{peri}$' + '\nperiplasmic density [fg/µm$^3$]', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/fig2_densities.pdf', bbox_inches='tight')