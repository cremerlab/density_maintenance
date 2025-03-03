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
fig, ax = plt.subplots(3, 2, figsize=(4, 3.1), sharex=True)
mapper = {'cyt_tot_per_cell': [ax[0, 0], 'black'],
          'rho_cyt_tot': [ax[0, 1], 'black'],
          'peri_prot_per_cell': [ax[2, 0], 'purple'],
          'rho_peri': [ax[2, 1], 'purple'],
          'mem_prot_per_cell': [ax[1, 0], 'blue'],
          'sigma_mem': [ax[1, 1], 'blue']}
# Set the range to compute the empirical trend
lam_range = np.linspace(0, 2.5)

# Iterate through each quantity and plot
res = {}
for g, d in densities[densities['quantity'].isin(mapper.keys())].groupby('quantity'):
    axis, color = mapper[g]
    fmt = size.viz.style_point('This Study')
    fmt['markeredgecolor'] = cor[f'primary_{color}']
    fmt['markerfacecolor'] = cor[f'pale_{color}']
    fmt['alpha'] = 0.75 
    axis.vlines(d['growth_rate_hr'], d['sig2_lower'], d['sig2_upper'], lw=0.5, 
                color=cor[f'primary_{color}'])
    axis.vlines(d['growth_rate_hr'], d['sig1_lower'], d['sig1_upper'], lw=1, 
                color=cor[f'primary_{color}'])
    axis.plot(d['growth_rate_hr'], d['mean'], **fmt)

    # Plot the empirical fits:
    if ('cyt_tot_per_cell' ==  g) or ('mem_prot_per_cell' == g) or (g == 'rho_peri'): 
        popt = scipy.stats.linregress(d['growth_rate_hr'], np.log(d['mean']))
        fit = np.exp(popt[1] + popt[0] * lam_range)
    else:
        popt = scipy.stats.linregress(d['growth_rate_hr'], d['mean'])
        fit = popt[1] + popt[0] * lam_range
    res[g] = popt 
    axis.plot(lam_range, fit, '--', color=cor[f'primary_{color}'], lw=1)
 
# Add context
for i in range(2):
    ax[-1, i].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0, 0].set_ylabel('$M_{prot}^{(cyt)} + M_{RNA}$\n[fg / cell]', fontsize=6)
ax[2, 0].set_ylabel('$M_{prot}^{(peri)}$\n[fg/ cell]', fontsize=6)
ax[1, 0].set_ylabel('$M_{prot}^{(mem)}$\n[fg/ cell]', fontsize=6)
ax[0, 1].set_ylabel(r'$\rho_{cyto}$' + '\n[fg / µm$^3$]', fontsize=6)
ax[2, 1].set_ylabel(r'$\rho_{peri}$' + '\n[fg / µm$^3$]', fontsize=6)
ax[1, 1].set_ylabel(r'$\sigma_{mem}$' + '\n[fg / µm$^2$]', fontsize=6)

# Control bounds
ax[0, 0].set_ylim([100, 1200])
ax[0, 1].set_ylim([0, 800])
ax[2, 0].set_ylim([0, 40])
ax[2, 1].set_ylim([30, 300])
ax[1, 0].set_ylim([0, 80])
ax[1, 1].set_ylim([0, 5])
plt.subplots_adjust(hspace=0.1)
plt.savefig('./plots/fig2_masses_densities.pdf', bbox_inches='tight')


#%% Plot the empirical kappa
fig, ax = plt.subplots(1, 1, figsize=(2, 1.3))
emp_kappa = densities[densities['quantity']=='empirical_kappa']
fmt = size.viz.style_point('This Study')
fmt['markeredgecolor'] = cor['primary_red']
fmt['markerfacecolor'] = cor['pale_red']
fmt['alpha'] = 0.75
ax.vlines(emp_kappa['growth_rate_hr'], emp_kappa['sig2_lower'], 
          emp_kappa['sig2_upper'], lw=0.5, color=cor['primary_red'])
ax.vlines(emp_kappa['growth_rate_hr'], emp_kappa['sig1_lower'], 
          emp_kappa['sig1_upper'], lw=1, color=cor['primary_red'])
ax.plot(emp_kappa['growth_rate_hr'], emp_kappa['mean'], **fmt)

# Compute and plot the empirical fit on the means
popt = scipy.stats.linregress(emp_kappa['growth_rate_hr'],
                              emp_kappa['mean'])
lam_range = np.linspace(0, 2.5)
fit = popt[1] + popt[0] * lam_range
ax.plot(lam_range, fit, '--', color=cor['primary_red'], lw=1)

# Add context
ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('cytoplasm-membrane\ndensity ratio [µm$^{-1}$]',fontsize=6)
ax.set_ylim([50, 250])
ax.set_xlim([0, 2.5])
plt.savefig('./plots/fig2_empirical_kappa.pdf', bbox_inches='tight')