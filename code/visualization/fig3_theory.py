#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the experimental data and restrict to wildtype
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data = data[data['strain']=='wildtype']

# Load the inference
preds = pd.read_csv('../../data/mcmc/theory_inference_ppc_summaries.csv')
pars = pd.read_csv('../../data/mcmc/theory_inference_parameter_summaries.csv')

# Load the empirical densities
densities = pd.read_csv('../../data/mcmc/empirical_densities_summary.csv')

#%% Plot the fit of phi_mem and phi_peri vs phi_rib
fig, ax = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)

# Plot our data
fmt = size.viz.style_point('This Study')
fmt['markeredgecolor'] = cor['primary_blue']
fmt['markerfacecolor'] = cor['pale_blue']
ax[0].plot(data['phi_rib'], data['phi_mem'], zorder=1000, **fmt)
fmt['markeredgecolor'] = cor['primary_purple']
fmt['markerfacecolor'] = cor['pale_purple']
ax[1].plot(data['phi_rib'], data['phi_peri'], zorder=1000, **fmt)

# Plot the PPCs
colors = ['blue', 'purple']
for i, q in enumerate(['phi_mem_ppc', 'phi_peri_ppc']):
    _pred = preds[preds['quantity']==q]
    ax[i].fill_between(_pred['phi_rib'], _pred['sig2_lower'], _pred['sig2_upper'],
                       color=cor[f'pale_{colors[i]}'], alpha=0.5)
    ax[i].fill_between(_pred['phi_rib'], _pred['sig1_lower'], _pred['sig1_upper'],
                       color=cor[f'light_{colors[i]}'], alpha=0.5) 
    ax[i].plot(_pred['phi_rib'], _pred['mean'], lw=1, color=cor[f'{colors[i]}'])
 
# Set context
ax[0].set_xlim([0.1, 0.4])
ax[0].set_ylim([0, 0.2])
ax[1].set_ylim([0, 0.15])
for a in ax:
    a.set_xlabel('ribosomal proteome allocation\n$\phi_{rib}$',
                 fontsize=6)
    a.set_xticks([0.1, 0.2, 0.3, 0.4])
ax[0].set_yticks([0, 0.05, 0.1, 0.15, 0.2])
ax[1].set_yticks([0, 0.05, 0.1, 0.15])

ax[0].set_ylabel('$\psi_{mem}$\nmembrane proteome parition', fontsize=6)
ax[1].set_ylabel('$\psi_{peri}$\nperiplasm proteome partition', fontsize=6)
plt.savefig('./plots/fig3_partition_trends.pdf', bbox_inches='tight')


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
ax.set_xlabel('ribosomal \nproteome allocation\n$\phi_{rib}$', fontsize=6)
ax.set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
plt.savefig('./plots/fig3_theory.pdf', bbox_inches='tight')