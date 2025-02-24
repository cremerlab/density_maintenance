#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
import scipy.stats
cor, pal = size.viz.matplotlib_style()

# Load the various datasets
rp_data = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')
rp_fit = pd.read_csv('../../data/mcmc/rna_protein_per_cell_ppc_summary.csv')
empirics = pd.read_csv('../../data/mcmc/empirical_densities_summary.csv')
empirics = empirics[empirics['strain']=='wildtype']

#%% Set figure canvas for RP fit
fig, ax = plt.subplots(2, 1, figsize=(2.5, 2.5), sharex=True)
fmt = size.viz.style_point('This Study')
fmt['marker'] = 'D'

# Plot the posterior predictive checks
for i, q in enumerate(['prot_per_cell_ppc', 'rna_per_cell_ppc']):
    d = rp_fit[rp_fit['quantity']==q]
    ax[i].fill_between(d['growth_rate_hr'], d['sig2_lower'], d['sig2_upper'], 
                       color=cor['pale_black'], alpha=0.75)
    ax[i].fill_between(d['growth_rate_hr'], d['sig1_lower'], d['sig1_upper'], 
                       color=cor['light_black'], alpha=0.75)
    ax[i].plot(d['growth_rate_hr'], d['mean'], '-', lw=1, color=cor['primary_black'])

# Plot our measurements
ax[0].plot(rp_data['growth_rate_hr'], rp_data['fg_protein_per_cell'], **fmt)
ax[1].plot(rp_data['growth_rate_hr'], rp_data['fg_rna_per_cell'], **fmt)

# Adjust scaling and add context
for a in ax:
    a.set_yscale('log')
ax[0].set_ylim([50, 1000])
ax[1].set_ylim([10, 1000])
ax[0].set_ylabel('$M_{prot}^{(tot)}$ [fg /cell]\ntotal protein per cell', fontsize=6)
ax[1].set_ylabel('$M_{RNA}^{(tot)}$ [fg /cell]\ntotal RNA per cell', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

plt.savefig('./plots/figS4_protein_rna_fits.pdf', bbox_inches='tight')

#%% Plot the compartment masses
fig, ax = plt.subplots(3, 1, figsize=(2.5, 3.3), sharex=True)

# Set colors and point styling
colors = ['black', 'purple', 'blue']
fmt = size.viz.style_point('This Study')

# Set the growth rate range for the empirical trend
lam_range = np.linspace(0, 2.5, 100)
for i, q in enumerate(['cyt_tot_per_cell', 'peri_prot_per_cell', 'mem_prot_per_cell']):
    d = empirics[empirics['quantity']==q]

    # Set the point context
    c = colors[i]
    fmt['markerfacecolor'] = cor[f'pale_{c}'] 
    fmt['markeredgecolor'] = cor[f'primary_{c}']

    # Plot the points
    ax[i].vlines(d['growth_rate_hr'], d['sig2_lower'], d['sig2_upper'], lw=0.5,
                 color=cor[f'primary_{c}'])
    ax[i].vlines(d['growth_rate_hr'], d['sig1_lower'], d['sig1_upper'], lw=1,
                 color=cor[f'primary_{c}'])
    ax[i].plot(d['growth_rate_hr'], d['mean'], **fmt)

    # Perform a regression
    if i != 1:
        popt = scipy.stats.linregress(d['growth_rate_hr'], np.log(d['mean']))
        fit = np.exp(popt[1] + popt[0] * lam_range)
    else:
        popt = scipy.stats.linregress(d['growth_rate_hr'], d['mean'])
        fit = popt[1] + popt[0] * lam_range
    ax[i].plot(lam_range, fit, '--', color=cor[f'primary_{c}'], lw=1)

# Set context
ax[0].set_ylim([0, 1000])
ax[0].set_ylabel('$\phi_{cyt}M_{prot}^{(tot)} + M_{RNA}^{(tot)}$ [fg / cell]\ncytoplasmic\nRNA + protein mass',
                 fontsize=6)
ax[1].set_ylim([0, 50])
ax[1].set_ylabel('$\phi_{peri}M_{prot}^{(tot)}$ [fg / cell]\nperiplasmic\nprotein mass', 
                 fontsize=6)
ax[2].set_ylim([0, 100])
ax[2].set_ylabel('$\phi_{mem}M_{prot}^{(tot)}$ [fg / cell]\nmembrane \nprotein mass',
                 fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/figS4_compartment_masses.pdf', bbox_inches='tight')