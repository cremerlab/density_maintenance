#%%
import numpy as np
import pandas as pd 
import scipy.stats
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

# Set constants
BETA = 1/0.4558
W_PERI = 0.025


# Load the data
data = pd.read_csv('../../data/compiled_measurements.csv')
data = data[data['strain']=='wildtype']

# Empirically calculate the density ratio for each point
data['emp_kappa'] = (2 * data['phi_cyto']/data['phi_mem']) / (data['sav_inv_um']**-1 - W_PERI)

# Load the parameter summarys
pars = pd.read_csv('../../data/mcmc/fig3_parameter_summaries.csv')

# Load the trend quantiles
quants = pd.read_csv('../../data/mcmc/fig3_fits.csv')

#%% Compute what's needed for the theory plot
# Set a range of kappas and conversion factor for phi_rib to R/P
kappa = [70, 90, 150, 175] 

# Perform regression on the trends to show contours of the theory
# Note that uncertainty bands exclusively come from Bayesian inferential model.
phi_mem_popt = scipy.stats.linregress(data['phi_rib'].values, data['phi_mem'].values)
phi_peri_popt = scipy.stats.linregress(data['phi_rib'], np.log(data['phi_peri'].values))

# Set a phi_rib range and compute the phi_mem and phi_peri
phi_rib_range = np.linspace(0, 0.45, 200)
phi_mem_range = phi_mem_popt[1] + phi_mem_popt[0] * phi_rib_range
phi_peri_range = np.exp(phi_peri_popt[1] + phi_peri_popt[1] * phi_rib_range)

# Compute the theory for the different kappas.
dfs = []
for kap in kappa:
    theory = kap * phi_mem_range / (2 * (1 + BETA * phi_rib_range - phi_mem_range - phi_peri_range)) 
    _df = pd.DataFrame(np.array([phi_rib_range, theory]).T, 
                       columns=['phi_rib', 'sav_theory'])
    _df['kappa'] = kap
    dfs.append(_df)
theory_df = pd.concat(dfs, sort=False)
theory_df

#%%
# Plot of the trends and the fits
fig, ax = plt.subplots(1, 2, figsize=(3, 1))

# Set postitional and stying information for the quantiles
axes = {'phi_mem_ppc': [ax[0], 'blue'],
        'phi_peri_ppc': [ax[1], 'purple']}
perc_colors = {'2sig':'pale_', '1sig':'light_'}

# Iterate through each quantity and plot
for g, d in quants[quants['quantity'].isin(axes)].groupby('quantity'):
    a = axes[g][0]
    a.plot(d['phi_rib'], d['mean_val'], '-', lw=1, color=cor[axes[g][1]])
    for k, v in perc_colors.items():
        c = cor[f'{v}{axes[g][1]}'] 
        a.fill_between(d['phi_rib'], d[f'{k}_lower'], d[f'{k}_upper'],
                       color=c, alpha=0.5)

# Plot markers
fmt = size.viz.style_point('This Study')
fmt['markeredgecolor'] = cor['primary_blue']
ax[0].plot(data['phi_rib'], data['phi_mem'], **fmt)
fmt['markeredgecolor'] = cor['primary_purple']
ax[1].plot(data['phi_rib'], data['phi_peri'], **fmt)

# Add context
for a in ax:
    a.set_xlabel('ribosomal proteome\nfraction\n$\phi_{rib}$', fontsize=6)
ax[0].set_ylabel('$\phi_{mem}$\nmembrane\nproteome fraction', fontsize=6) 
ax[1].set_ylabel('$\phi_{peri}$\nperiplasmic\nproteome fraction', fontsize=6)

# Adjust axes
ax[0].set_xlim([0.1, 0.3])
ax[1].set_xlim([0.1, 0.3])
ax[0].set_ylim([0, 0.20])
ax[1].set_ylim([0, 0.15])

plt.savefig('./plots/fig3.2_trends.pdf', bbox_inches='tight')

#%% Plot the theory-experiment agreement
fig, ax = plt.subplots(1, 1, figsize=(2,2))

# Plot the lines for kappa and add labels afterward
for g, d in theory_df.groupby('kappa'):
    ax.plot(d['phi_rib'], d['sav_theory'], '-', lw=0.75, color=cor['pale_black'])

# Plot the fit and credible regions
fit = quants[quants['quantity']=='theory_ppc']
ax.fill_between(fit['phi_rib'], fit['2sig_lower'], fit['2sig_upper'], color=cor['pale_green'], alpha=0.75)
ax.fill_between(fit['phi_rib'], fit['1sig_lower'], fit['1sig_upper'], color=cor['light_green'])
ax.plot(fit['phi_rib'], fit['mean_val'],'-', lw=1, color=cor['green'])

# Plot our experiments
fmt = size.viz.style_point('This Study')
ax.plot(data['phi_rib'], data['sav_inv_um'], **fmt)

# Add context
ax.set_xlim([0.08, 0.35])
ax.set_ylim([3, 10])
ax.set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax.set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
plt.savefig('./plots/fig3.2_sav_theory.pdf', bbox_inches='tight')

#%% Plot the agreement with the empirically calculated kappa 
fig, ax = plt.subplots(1, 1, figsize=(1.5, 0.75))

# Plot the credible regions
kappa_est = pars[pars['quantity']=='kappa']
ax.fill_between([0, 3], kappa_est['2sig_lower'], kappa_est['2sig_upper'], 
                color=cor['pale_red'], alpha=0.5)
ax.fill_between([0, 3], kappa_est['1sig_lower'], kappa_est['1sig_upper'], 
                color=cor['light_red'])
ax.hlines(kappa_est['mean_val'], 0, 3, lw=1, color=cor['red'])

# Plot the markers
fmt = size.viz.style_point('This Study')
ax.plot(data['growth_rate_hr'], data['emp_kappa'], **fmt)


 
# Add context
ax.set_ylim([50, 180])
ax.set_xlim([0.1, 2.5])
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
