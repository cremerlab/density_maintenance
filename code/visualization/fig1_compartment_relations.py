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
ax[0].plot(data['growth_rate_hr'], data['psi_cyto'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['psi_cyto'])
print(f'psi_cyto slope: {popt[0]:0.3f} ({popt[0] /data["psi_cyto"].mean():0.3f})')
print(f'psi_cyto intercept: {popt[1]:0.3f} ({popt[1] /data["psi_cyto"].mean():0.3f})')
ax[0].plot(lam_range, popt[1] + popt[0] * lam_range, '--',
           color=cor['primary_black'], lw=1)

# Plot the periplasmic allocation
fmt['markeredgecolor'] = cor['primary_purple']
fmt['markerfacecolor'] = cor['pale_purple']
ax[1].plot(data['growth_rate_hr'], data['psi_peri'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['psi_peri'])
print(f'psi_peri slope: {popt[0]:0.3f} ({popt[0] / data["psi_peri"].mean():0.3f})')
print(f'psi_peri intercept: {popt[1]:0.3f} ({popt[1] / data["psi_peri"].mean():0.3f})')
ax[1].plot(lam_range, popt[1] + popt[0] * lam_range, '--',
           color=cor['primary_purple'], lw=1)

# Plot membrane allocation
fmt['markerfacecolor'] = cor['pale_blue']
fmt['markeredgecolor'] = cor['primary_blue']
ax[2].plot(data['growth_rate_hr'], data['psi_mem'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['psi_mem'])
print(f'psi_mem slope: {popt[0]:0.3f} ({popt[0] / data["psi_mem"].mean():0.3f})')
print(f'psi_mem intercept: {popt[1]:0.3f} ({popt[1] / data["psi_mem"].mean():0.3f})')
ax[2].plot(lam_range, popt[1] + popt[0] * lam_range, '--',
           color=cor['primary_blue'], lw=1)


# Add context
ax[0].set_ylabel('$\psi_{cyto}$\ncytoplasmic\nproteome partition\n', fontsize=6)
ax[1].set_ylabel('$\psi_{peri}$\nperiplasm\nproteome partition\n', fontsize=6)
ax[2].set_ylabel('$\psi_{mem}$\nmembrane\nproteome partition\n', fontsize=6)
ax[2].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

# Set axis limits
ax[0].set_ylim([0.7, 0.9])
ax[0].set_yticks([0.7, 0.75, 0.8, 0.85, 0.9])
ax[1].set_ylim([0, 0.15])
ax[1].set_yticks([0, 0.05, 0.1, 0.15])
ax[2].set_ylim([0, 0.2])
ax[2].set_yticks([0, 0.05, 0.1, 0.15, 0.2])
plt.savefig('./plots/fig1_compartment_partition.pdf', bbox_inches='tight')

#%%  Plot (E - F)
fig, ax = plt.subplots(2, 1, figsize=(1.6, 3))
fmt = size.viz.style_point('This Study')

# Plot the tradoff between periplasm and cytoplasm
ax[0].plot(data['psi_peri'], data['psi_cyto'], **fmt)
# Plot a line with a slope of -1 showing a strong tradeoff
psi_peri_range = np.array([0.05, 0.11])
psi_cyto_range = np.array([0.85, 0.79])
ax[0].plot(psi_peri_range, psi_cyto_range, '-', color=cor['light_black'], 
           lw=0.5)

# plot the conservation of mass between compartments
ax[1].plot(data['growth_rate_hr'], data['psi_peri'] + data['psi_cyto'], **fmt)
popt = scipy.stats.linregress(data['growth_rate_hr'], data['psi_peri'].values + data['psi_cyto'].values)
print(f'psi_peri + psi_cyto slope: {popt[0]:0.3f} ({popt[0] / (data["psi_peri"] + data["psi_cyto"]).mean():0.3f})')
print(f'psi_peri + psi_cyto intercept: {popt[1]:0.3f} ({popt[1] / (data["psi_peri"] + data["psi_cyto"]).mean():0.3f})')
ax[1].plot(lam_range, popt[1] + popt[0] * lam_range, 'k--', lw=1)

# Set context
ax[0].set_xlabel('periplasm proteome partition\n$\psi_{peri}$', fontsize=6)
ax[0].set_ylabel('$\psi_{cyto}$\ncytoplasm\nproteome partition', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[1].set_ylabel('$\psi_{cyto} + \psi_{peri}$\nperiplasm + cytoplasm\nproteome partition', fontsize=6)

# Set limits and adjust subplots
ax[0].set_ylim([0.75, 0.9])
ax[0].set_xlim([0, 0.15])
ax[1].set_ylim([0.6, 1.0])
plt.subplots_adjust(hspace=0.35)
plt.savefig('./plots/fig1_cytoplasm_periplasm_tradeoff.pdf', bbox_inches='tight')
