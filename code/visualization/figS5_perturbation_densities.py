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
fig, ax = plt.subplots(2,2, figsize=(3, 3))
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
ax[0].plot(size_data['growth_rate_hr'], size_data['volume_fL'], **fmt)

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

#%%
fig, ax = plt.subplots(1,1, figsize=(2, 1))

perturbs['calc_kappa'] = ((perturbs['phi_cyto'] + perturbs['phi_rib'] * (1/0.4558)) / perturbs['phi_mem']) * ((2 * perturbs['surface_area_um2'])/(perturbs['volume_fL'] - 0.025 * perturbs['surface_area_um2']))
wt_data['calc_kappa'] = ((wt_data['phi_cyto'] + wt_data['phi_rib'] * (1/0.4558)) / wt_data['phi_mem']) * ((2 * wt_data['surface_area_um2'])/(wt_data['volume_fL'] - 0.025 * wt_data['surface_area_um2']))

#
rel = perturbs[(perturbs['strain']=='relA') & (perturbs['carbon_source']=='glucose')]
ax.plot(wt_data['growth_rate_hr'], wt_data['calc_kappa'], 'o', color=cor['pale_black'], ms=4)
ax.plot(perturbs['growth_rate_hr'], perturbs['calc_kappa'], 'o', ms=4)
ax.plot(rel['growth_rate_hr'], rel['calc_kappa'], 'rD', ms=4)

ax.set_ylim([50, 200])

