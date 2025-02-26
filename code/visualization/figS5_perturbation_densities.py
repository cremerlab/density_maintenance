#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the various data sets
prot_data = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')
size_data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
wt_data = size_data[size_data['strain']=='wildtype']
perturbs = size_data[size_data['strain']!='wildtype']
wt_kappas = pd.read_csv('../../data/mcmc/empirical_densities_summary.csv')
perturb_kappas = pd.read_csv('../../data/mcmc/perturbation_empirical_densities_summary.csv')
density = pd.read_csv('../../data/mcmc/estimated_constant_total_protein_density.csv')
ppcs = pd.read_csv('../../data/mcmc/volume_protein_per_cell_ppc_summary.csv')

# Set the figure canvas for the empirical fits.
fig, ax = plt.subplots(2,2, figsize=(3.25, 3))
ax[0, 0].set_visible(False)
ax = [ax[0,1], ax[1, 0], ax[1, 1]]

# Format point
fmt = size.viz.style_point('This Study')

# Plot volume ppc
vol_ppc = ppcs[ppcs['quantity']=='volume_ppc']
ax[0].fill_between(vol_ppc['growth_rate_hr'], vol_ppc['sig2_lower'], 
                   vol_ppc['sig2_upper'], color=cor['pale_black'], alpha=0.75)
ax[0].fill_between(vol_ppc['growth_rate_hr'], vol_ppc['sig1_lower'], 
                   vol_ppc['sig1_upper'], color=cor['pale_black'], alpha=0.75)
ax[0].plot(vol_ppc['growth_rate_hr'], vol_ppc['mean'], '-', lw=1, color=cor['primary_black'])
ax[0].plot(wt_data['growth_rate_hr'], wt_data['volume_fL'], **fmt)

# Plot distribution of rho_prot_mu
ax[1].hist(density['rho_prot_mu'], bins=35, color=cor['light_black'],
           edgecolor=cor['primary_black'], linewidth=0.5)


# Plot the protein per cell data and the empirical fits
fmt = size.viz.style_point('This Study')
fmt['marker'] = 'D'
prot_ppc = ppcs[ppcs['quantity']=='protein_per_cell_ppc']
ax[2].fill_between(prot_ppc['growth_rate_hr'], prot_ppc['sig2_lower'], 
                   prot_ppc['sig2_upper'], color=cor['pale_black'], alpha=0.75)
ax[2].fill_between(prot_ppc['growth_rate_hr'], prot_ppc['sig1_lower'], 
                   prot_ppc['sig1_upper'], color=cor['pale_black'], alpha=0.75)
ax[2].plot(prot_ppc['growth_rate_hr'], prot_ppc['mean'], '-', lw=1, color=cor['primary_black'])
ax[2].plot(prot_data['growth_rate_hr'], prot_data['fg_protein_per_cell'], **fmt)


# Add context
ax[0].set_yscale('log')
ax[0].set_ylim([0.1, 10])
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$V$\ncell volume [µm$^3$]', fontsize=6)

ax[1].set_yticks([])
ax[1].set_xlabel('total protein density [fg/µm$^3$]\n' + r'$\rho_{prot}^{(tot)}$',
                 fontsize=6)
ax[2].set_yscale('log')
ax[2].set_ylabel('$M_{prot}^{(tot)}$\n total protein [fg / cell]', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/figS5_inference.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(4,1, figsize=(2.5, 4), sharex=True)

# Set axes and styles
axes = {'rho_cyt_tot': [ax[0], 'black'], 
        'rho_peri': [ax[1], 'purple'], 
        'sigma_mem': [ax[2], 'blue'],
        'empirical_kappa': [ax[3], 'red']}
fmt = size.viz.style_point('This Study')
fmt['markeredgewidth'] = 0
fmt['alpha'] = 0.4

# Plot the wildtype data
for q, (a, c) in axes.items():
    wt = wt_kappas[wt_kappas['quantity']==q] 
    fmt['markerfacecolor'] = cor[f'pale_black']
    a.vlines(wt['growth_rate_hr'], wt['sig2_lower'], wt['sig2_upper'], lw=0.5, 
             color=cor[f'pale_black'], alpha=fmt['alpha'])
    a.vlines(wt['growth_rate_hr'], wt['sig1_lower'], wt['sig1_upper'], lw=1, 
             color=cor[f'pale_black'], alpha=fmt['alpha'])
    a.plot(wt['growth_rate_hr'], wt['mean'], **fmt)

# Plot the rel perturbations
rel = perturb_kappas[perturb_kappas['strain']=='relA']
colors = {0:'pale_', 2:'primary_', 4:''}
markers = {'glucose': 'D', 'glucoseCAA': 's'}
fmt = size.viz.style_point('This Study')
fmt['markersize'] = 3
fmt['alpha'] = 0.75
axes = {'rho_cyt': ax[0], 
        'rho_peri': ax[1], 
        'sigma_mem': ax[2],
        'empirical_kappa': ax[3]}
for g, d in rel.groupby(['carbon_source', 'inducer_conc']):
    for q, a in axes.items():
        _d = d[d['quantity']==q]
        fmt['markerfacecolor'] = cor[f'{colors[g[1]]}red']
        fmt['markeredgecolor'] = cor['dark_red']
        fmt['marker'] = markers[g[0]]
        if q != 'empirical_kappa':  # otherwise, directly calculated with no uncertainty.   
            a.vlines(_d['growth_rate_hr'], _d['sig2_lower'], _d['sig2_upper'], lw=0.5,
                 color=cor[f'dark_red'], alpha=fmt['alpha'])
            a.vlines(_d['growth_rate_hr'], _d['sig1_lower'], _d['sig1_upper'], lw=1,
                 color=cor[f'dark_red'], alpha=fmt['alpha'])
        a.plot(_d['growth_rate_hr'], _d['mean'], **fmt)

mesh = perturb_kappas[perturb_kappas['strain']=='meshI']
colors = {0:'pale_', 100:''}
markers = {'glucose': 'D', 'glucoseCAA': 's'}
fmt = size.viz.style_point('This Study')
fmt['markersize'] = 3
fmt['alpha'] = 0.75
axes = {'rho_cyt': ax[0], 
        'rho_peri': ax[1], 
        'sigma_mem': ax[2],
        'empirical_kappa': ax[3]}
for g, d in mesh.groupby(['carbon_source', 'inducer_conc']):
    for q, a in axes.items():
        _d = d[d['quantity']==q]
        fmt['markerfacecolor'] = cor[f'{colors[g[1]]}gold']
        fmt['markeredgecolor'] = cor['dark_gold']
        fmt['marker'] = markers[g[0]]
        if q != 'empirical_kappa':  # otherwise, directly calculated with no uncertainty.
            a.vlines(_d['growth_rate_hr'], _d['sig2_lower'], _d['sig2_upper'], lw=0.5,
                 color=cor[f'dark_gold'], alpha=fmt['alpha'])
            a.vlines(_d['growth_rate_hr'], _d['sig1_lower'], _d['sig1_upper'], lw=1,
                 color=cor[f'dark_gold'], alpha=fmt['alpha'])
        a.plot(_d['growth_rate_hr'], _d['mean'], **fmt)

# Add context
ax[0].set_ylim([100, 600])
ax[0].set_ylabel(r'$\rho_{cyt}$' + '\ncytoplasmic\ndensity [fg/µm$^3$]', fontsize=6)

ax[1].set_ylim([0, 300])
ax[1].set_ylabel(r'$\rho_{peri}$' + '\nperiplasmic \ndensity [fg/µm$^3$]', fontsize=6)

ax[2].set_ylim([0, 5])
ax[2].set_ylabel('$\sigma_{mem}$\n membrane\ndensity [fg / µm$^2$]', fontsize=6)

ax[3].set_ylim([50, 250])
ax[3].set_xlim([0.25, 1.5])
ax[3].set_ylabel('cytoplasm-membrane\ndensity ratio [µm$^{-1}$]', fontsize=6)
ax[3].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
plt.savefig('./plots/figS5_densities.pdf', bbox_inches='tight')