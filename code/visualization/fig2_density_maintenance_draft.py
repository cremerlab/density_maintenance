#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
import scipy.stats
cor, pal = size.viz.matplotlib_style()

# Load the experimental data and restrict to wildtype
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data = data[data['strain']=='wildtype']

# Load the inference
preds = pd.read_csv('../../data/mcmc/theory_inference_ppc_summaries.csv')
pars = pd.read_csv('../../data/mcmc/theory_inference_parameter_summaries.csv')

# Load the empirical densities
densities = pd.read_csv('../../data/mcmc/empirical_densities_summary.csv')

#%% Plot the masses and densities within the compartments
fig, ax = plt.subplots(3, 2, figsize=(5.5, 3.1), sharex=True)
mapper = {'cyt_tot_per_cell': [ax[0, 0], cor['primary_black']],
          'rho_cyt_tot': [ax[0, 1], cor['primary_black']],
          'peri_prot_per_cell': [ax[1, 0], cor['primary_purple']],
          'rho_peri': [ax[1, 1], cor['primary_purple']],
          'mem_prot_per_cell': [ax[2, 0], cor['primary_blue']],
          'sigma_mem': [ax[2, 1], cor['primary_blue']]}
# Set the range to compute the empirical trend
lam_range = np.linspace(0, 2.5)

# Iterate through each quantity and plot
res = {}
for g, d in densities[densities['quantity'].isin(mapper.keys())].groupby('quantity'):
    axis, color = mapper[g]
    fmt = size.viz.style_point('This Study')
    fmt['markeredgecolor'] = color
    fmt['alpha'] = 0.75 
    axis.vlines(d['growth_rate_hr'], d['sig2_lower'], d['sig2_upper'], lw=0.5, 
                color=color)
    axis.vlines(d['growth_rate_hr'], d['sig1_lower'], d['sig1_upper'], lw=1, 
                color=color)
    axis.plot(d['growth_rate_hr'], d['mean'], **fmt)

    # Plot the empirical fits:
    if ('cyt_tot_per_cell' ==  g) or ('mem_prot_per_cell' == g) or (g == 'rho_peri'): 
        popt = scipy.stats.linregress(d['growth_rate_hr'], np.log(d['mean']))
        fit = np.exp(popt[1] + popt[0] * lam_range)
    else:
        popt = scipy.stats.linregress(d['growth_rate_hr'], d['mean'])
        fit = popt[1] + popt[0] * lam_range
    res[g] = popt 
    axis.plot(lam_range, fit, '--', color=color, lw=1)
 
# Add context
for i in range(2):
    ax[-1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0, 0].set_ylabel('$M_{prot}^{(cyt)} + M_{RNA}$\n[fg / cell]', fontsize=6)
ax[1, 0].set_ylabel('$M_{prot}^{(peri)}$\n[fg/ cell]', fontsize=6)
ax[2, 0].set_ylabel('$M_{prot}^{(mem)}$\n[fg/ cell]', fontsize=6)
ax[0, 1].set_ylabel(r'$\rho_{cyt}$' + '\n[fg / µm$^3$]', fontsize=6)
ax[1, 1].set_ylabel(r'$\rho_{peri}$' + '\n[fg / µm$^3$]', fontsize=6)
ax[2, 1].set_ylabel(r'$\sigma_{mem}$' + '\n[fg / µm$^2$]', fontsize=6)

# Control bounds
ax[0, 0].set_ylim([100, 1200])
ax[0, 1].set_ylim([100, 700])
ax[1, 0].set_ylim([0, 40])
ax[1, 1].set_ylim([30, 300])
ax[2, 0].set_ylim([0, 80])
ax[2, 1].set_ylim([0, 5])
plt.subplots_adjust(hspace=0.1)
plt.savefig('./plots/fig2_masses_densities.pdf', bbox_inches='tight')

#%% Plot the fit of phi_mem and phi_peri vs phi_rib
fig, ax = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)

# Plot our data
fmt = size.viz.style_point('This Study')
fmt['markeredgecolor'] = cor['primary_blue']
ax[0].plot(data['phi_rib'], data['phi_mem'], zorder=1000, **fmt)
fmt['markeredgecolor'] = cor['primary_purple']
ax[1].plot(data['phi_rib'], data['phi_peri'], zorder=1000, **fmt)

# Plot the PPCs
colors = ['blue', 'purple']
for i, q in enumerate(['phi_mem_ppc', 'phi_peri_ppc']):
    _pred = preds[preds['quantity']==q]
    ax[i].fill_between(_pred['phi_rib'], _pred['sig2_lower'], _pred['sig2_upper'],
                       color=cor[f'pale_{colors[i]}'], alpha=0.75)
    ax[i].fill_between(_pred['phi_rib'], _pred['sig1_lower'], _pred['sig1_upper'],
                       color=cor[f'light_{colors[i]}'], alpha=0.75) 
    ax[i].plot(_pred['phi_rib'], _pred['mean'], lw=1, color=cor[f'{colors[i]}'])
 
# Set context
ax[0].set_xlim([0.1, 0.4])
ax[0].set_ylim([0, 0.2])
ax[1].set_ylim([0, 0.15])
for a in ax:
    a.set_xlabel('ribosomal allocation\n$\phi_{rib}$',
                 fontsize=6)
    a.set_xticks([0.1, 0.2, 0.3, 0.4])
ax[0].set_yticks([0, 0.05, 0.1, 0.15, 0.2])
ax[1].set_yticks([0, 0.05, 0.1, 0.15])

ax[0].set_ylabel('$\phi_{mem}$\nmembrane protein allocation', fontsize=6)
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic protein allocation', fontsize=6)
plt.savefig('./plots/fig2_allocation_trends.pdf', bbox_inches='tight')


#%% Plot the SAV theory and estimated values of kappa. 
# Set a helper function to compute the theory for different kappa values. 
phi_rib_range = np.linspace(0.08, 0.4, 100)
def prediction_helper(kappa: np.ndarray) -> np.ndarray:
    # Extract the empirical parameters for fits
    beta_0_phi_mem = pars[pars['quantity']=='beta_0_phi_mem']['mean'].values[0]
    beta_1_phi_mem = pars[pars['quantity']=='beta_1_phi_mem']['mean'].values[0]
    beta_0_phi_peri = pars[pars['quantity']=='beta_0_phi_peri']['mean'].values[0]
    beta_1_phi_peri = pars[pars['quantity']=='beta_1_phi_peri']['mean'].values[0]
    
    # Compute the empirical allocation 
    phi_mem = beta_0_phi_mem + beta_1_phi_mem * phi_rib_range
    phi_peri = beta_0_phi_peri * np.exp(beta_1_phi_peri * phi_rib_range)

    # Compute the theory
    numer = kappa * phi_mem 
    denom = 2 * (1 + (1/0.4558) * phi_rib_range - phi_mem - phi_peri)
    return numer / denom


# Set the range of kappas 
kappas = [70, 90, 175, 200]

# Define the figure canvas 
fig, ax = plt.subplots(1, 1, figsize=(2, 2.5))
for k in kappas:
    ax.plot(phi_rib_range, prediction_helper(k), '-', lw=1, color=cor['pale_black'])

# Plot the inference
theo = preds[preds['quantity'] == 'theory_ppc']
ax.fill_between(theo['phi_rib'], theo['sig2_lower'], theo['sig2_upper'],
                color=cor['pale_green'], alpha=0.75)
ax.fill_between(theo['phi_rib'], theo['sig1_lower'], theo['sig1_upper'],
                color=cor['light_green'], alpha=0.75)
ax.plot(theo['phi_rib'], theo['mean'], lw=1, color=cor['primary_green'])

# Plot the wildtype data
fmt = size.viz.style_point('This Study')
ax.plot(data['phi_rib'], data['surface_to_volume_inv_um'], **fmt)

# Set context
ax.set_xlim([0.08, 0.35])
ax.set_ylim([3, 10])
ax.set_xlabel('ribosomal allocation\n$\phi_{rib}$', fontsize=6)
ax.set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
plt.savefig('./plots/fig2_theory.pdf', bbox_inches='tight')