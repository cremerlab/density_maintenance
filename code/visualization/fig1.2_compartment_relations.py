#%%
import numpy as np 
import pandas as pd 
import scipy.stats
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the experimental data and restrict to wildtype
data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data = data[data['strain']=='wildtype']

#%% Plot (B - D) the allocation across the three  compartments
fig, ax = plt.subplots(3, 1, figsize=(2.3, 3), sharex=True)
# Set the growth rate ratge over whidh to plot the empirical fits
lam_range = np.linspace(0, 2.5, 30)

# Plot cytoplasmic allocation
fmt = size.viz.style_point('This Study')
fmt['markerfacecolor'] = cor['pale_black']
fmt['alpha'] = 0.8
ax[0].plot(data['growth_rate_hr'], data['phi_cyto'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['phi_cyto'])
ax[0].plot(lam_range, popt[1] + popt[0] * lam_range, '--',
           color=cor['primary_black'], lw=1)

# Plot membrane allocation
fmt['markerfacecolor'] = cor['pale_blue']
fmt['markeredgecolor'] = cor['primary_blue']
ax[1].plot(data['growth_rate_hr'], data['phi_mem'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['phi_mem'])
ax[1].plot(lam_range, popt[1] + popt[0] * lam_range, '--',
           color=cor['primary_blue'], lw=1)

# Plot the periplasmic allocation
fmt['markeredgecolor'] = cor['primary_purple']
fmt['markerfacecolor'] = cor['pale_purple']
ax[2].plot(data['growth_rate_hr'], data['phi_peri'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], np.log(data['phi_peri']))
ax[2].plot(lam_range, np.exp(popt[1] + popt[0] * lam_range), '--',
           color=cor['primary_purple'], lw=1)

# Add context
ax[0].set_ylabel('$\phi_{cyto}$\ncytoplasmic\nproteome allocation\n', fontsize=6)
ax[1].set_ylabel('$\phi_{mem}$\nmembrane\nproteome allocation\n', fontsize=6)
ax[2].set_ylabel('$\phi_{peri}$\nperiplasm\nproteome allocation\n', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

# Set axis limits
ax[0].set_ylim([0.6, 1.0])
ax[0].set_yticks([0.60, 0.7, 0.8, 0.9, 1])
ax[0].set_yticklabels(['0.60', '0.70', '0.80', '0.90', '1.00'])
ax[1].set_ylim([0, 0.2])
ax[1].set_yticks([0, 0.05, 0.1, 0.15, 0.2])
ax[2].set_ylim([0, 0.15])
ax[2].set_yticks([0, 0.05, 0.1, 0.15])
plt.savefig('./plots/fig1.2_compartment_allocation.pdf', bbox_inches='tight')

#%%  Plot (E - F)
fig, ax = plt.subplots(2, 1, figsize=(1.6, 3))
fmt = size.viz.style_point('This Study')

# Plot the tradoff between periplasm and cytoplasm
ax[0].plot(data['phi_peri'], data['phi_cyto'], **fmt)
phi_peri_range = np.linspace(0, 0.2, 30)
popt = scipy.stats.linregress(data['phi_peri'], data['phi_cyto'])
ax[0].plot(phi_peri_range, popt[1] + popt[0] * phi_peri_range, 'k--', lw=1)

# plot the conservation of mass between compartments
ax[1].plot(data['growth_rate_hr'], data['phi_peri'] + data['phi_cyto'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['phi_peri'].values + data['phi_cyto'].values)
ax[1].plot(lam_range, popt[1] + popt[0] * lam_range, 'k--', lw=1)

# Set context
ax[0].set_xlabel('periplasm proteome allocation\n$\phi_{peri}$', fontsize=6)
ax[0].set_ylabel('$\phi_{cyto}$\ncytoplasm\nproteome allocation', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1].set_ylabel('$\phi_{cyto} + \phi_{peri}$\nperiplasm + cytoplasm\nproteome allocation', fontsize=6)

# Set limits and adjust subplots
ax[0].set_ylim([0.75, 0.9])
ax[0].set_xlim([0, 0.15])
ax[1].set_ylim([0.6, 1.0])
plt.subplots_adjust(hspace=0.35)
plt.savefig('./plots/fig1.2_cytoplasm_periplasm_tradeoff.pdf', bbox_inches='tight')
