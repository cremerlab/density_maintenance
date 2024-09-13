#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import size.viz 
np.random.seed(3049344) # Randomly set for the posterior predictive check
cor, pal = size.viz.matplotlib_style()

BETA = 1/0.4558

# Load our mass spec data
data = pd.read_csv('../processing/mass_spectrometry/total_collated_data.csv')
data = data[data['strain']!='lacZ']
wt = data[data['strain']=='wildtype']

# Load the results of the full inference
samples = pd.read_csv('../analysis/output/phi_rib_scaling_fits_summary.csv')
samples = samples[samples['interval'].isin(['95%', '68%', 'median'])]
kappa = pd.read_csv('../analysis/output/kappa_density_samples.csv')
phi_rib_samples = pd.read_csv('../analysis/output/phi_rib_scaling_fits_samples.csv')

# Load the literature mass spec data
lit_data = pd.read_csv('../../data/collated/literature_mass_spec_aggregated.csv')

# Set a mapper for colors
mapper = {'95%':'pale_', '68%':'primary_', 'median':''}

fig, ax = plt.subplots(2, 2, figsize=(4, 3))
ax = ax.ravel()
ax[0].set_ylim([5, 8])
ax[0].set_xlim([0.1, 0.2])

# Plot the prediction
pred = samples[samples['quantity']=='sav_ppc']
for g, d in pred.groupby('interval', sort=False):
    # Do the fill between
    if g == 'median':
        ax[0].plot(d['phi_rib'], d['lower'], color=cor[f'{mapper[g]}green'], lw=1)
    else:
        ax[0].fill_between(d['phi_rib'], d['lower'], d['upper'], color=cor[f'{mapper[g]}green'], alpha=0.4)

# plot the wildtype samples
fmt = size.viz.style_point('This Study')
fmt['alpha'] = 0.5
ax[0].plot(wt['phi_rib'], wt['surface_to_volume'], **fmt)

# Compute the predicted versus measured SAV
fmt = size.viz.style_point('This Study')
fmt['markersize'] = 3
fmt['alpha'] = 0.5
color_mapper = {0:'light_',  2: 'primary_', 4:'', 100:''}

# Plot the scaling predicted vs measured
for i in range(len(wt)):
    pred_phi_mem_mu = phi_rib_samples.beta_0_phi_mem + phi_rib_samples.beta_1_phi_mem * wt['phi_rib'].values[i]
    pred_phi_mem = np.random.normal(pred_phi_mem_mu, phi_rib_samples.phi_mem_sigma)
    pred_median = np.median(pred_phi_mem)
    pred_2sig =  np.percentile(pred_phi_mem, [2.5, 97.5])
    pred_sig = np.percentile(pred_phi_mem, [16, 84]) 
    ax[1].vlines(wt['phi_mem'].values[i], pred_2sig[0], pred_2sig[1], lw=0.5, color=fmt['markeredgecolor'])
    ax[1].vlines(wt['phi_mem'].values[i], pred_sig[0], pred_sig[1], lw=1, color=fmt['markeredgecolor'])
    ax[1].plot(wt['phi_mem'].values[i], pred_median, **fmt)

    pred_phi_peri_mu = phi_rib_samples.beta_0_phi_peri + phi_rib_samples.beta_1_phi_peri * wt['phi_rib'].values[i]
    pred_phi_peri = np.exp(np.random.normal(pred_phi_peri_mu, phi_rib_samples.phi_peri_sigma))
    pred_median = np.median(pred_phi_peri)
    pred_2sig =  np.percentile(pred_phi_peri, [2.5, 97.5])
    pred_sig = np.percentile(pred_phi_peri, [16, 84]) 
    ax[2].vlines(wt['phi_peri'].values[i], pred_2sig[0], pred_2sig[1], lw=0.5, color=fmt['markeredgecolor'])
    ax[2].vlines(wt['phi_peri'].values[i], pred_sig[0], pred_sig[1], lw=1, color=fmt['markeredgecolor'])
    ax[2].plot(wt['phi_peri'].values[i], pred_median, **fmt) 

ax[1].plot([0.05, 0.18], [0.05, 0.18], '--', color=cor['primary_black'])
ax[1].set_ylim([0.05, 0.18])



for i in range(len(wt)):
   pred_SAV_mu = kappa.kappa * wt['phi_mem'].values[i] / (2 * (1 + BETA * wt['phi_rib'].values[i] - wt['phi_mem'].values[i] - wt['phi_peri'].values[i]))
   pred_SAV = np.random.normal(pred_SAV_mu, kappa.sav_sigma)
   pred_median = np.median(pred_SAV)
   pred_2sig =  np.percentile(pred_SAV, [2.5, 97.5])
   pred_sig = np.percentile(pred_SAV, [16, 84])
   ax[3].vlines(wt['surface_to_volume'].values[i], pred_2sig[0], pred_2sig[1], lw=0.5, color=fmt['markeredgecolor'])
   ax[3].vlines(wt['surface_to_volume'].values[i], pred_sig[0], pred_sig[1], lw=1, color=fmt['markeredgecolor'])
   ax[3].plot(wt['surface_to_volume'].values[i], pred_median, **fmt)

for g, d in data[data['strain']!='wildtype'].groupby(['strain', 'inducer_conc']):
    if g[0] == 'relA':
       c = 'red' 
    else:
        c = 'gold'
    fmt = size.viz.style_point('This Study')
    fmt['markeredgecolor'] = cor[f'dark_{c}']
    fmt['color'] = cor[f'{color_mapper[g[1]]}{c}']
    ax[0].plot(d['phi_rib'], d['surface_to_volume'], **fmt)

    # Compute the predicted versus measured SAV
    for i in range(len(d)):
        pred_SAV_mu = kappa.kappa * d['phi_mem'].values[i] / (2 * (1 + BETA * d['phi_rib'].values[i] - d['phi_mem'].values[i] - d['phi_peri'].values[i]))
        pred_SAV = np.random.normal(pred_SAV_mu, kappa.sav_sigma)
        pred_median = np.median(pred_SAV)
        pred_2sig =  np.percentile(pred_SAV, [2.5, 97.5])
        pred_sig = np.percentile(pred_SAV, [16, 84])
        ax[3].vlines(d['surface_to_volume'].values[i], pred_2sig[0], pred_2sig[1], lw=0.5, color=fmt['markeredgecolor'])
        ax[3].vlines(d['surface_to_volume'].values[i], pred_sig[0], pred_sig[1], lw=1, color=fmt['markeredgecolor'])
        ax[3].plot(d['surface_to_volume'].values[i], pred_median, **fmt)

        pred_phi_mem_mu = phi_rib_samples.beta_0_phi_mem + phi_rib_samples.beta_1_phi_mem * d['phi_rib'].values[i]
        pred_phi_mem = np.random.normal(pred_phi_mem_mu, phi_rib_samples.phi_mem_sigma)
        pred_median = np.median(pred_phi_mem)
        pred_2sig =  np.percentile(pred_phi_mem, [2.5, 97.5])
        pred_sig = np.percentile(pred_phi_mem, [16, 84]) 
        ax[1].vlines(d['phi_mem'].values[i], pred_2sig[0], pred_2sig[1], lw=0.5, color=fmt['markeredgecolor'])
        ax[1].vlines(d['phi_mem'].values[i], pred_sig[0], pred_sig[1], lw=1, color=fmt['markeredgecolor'])
        ax[1].plot(d['phi_mem'].values[i], pred_median, **fmt) 


        pred_phi_peri_mu = phi_rib_samples.beta_0_phi_peri + phi_rib_samples.beta_1_phi_peri * d['phi_rib'].values[i]
        pred_phi_peri = np.exp(np.random.normal(pred_phi_peri_mu, phi_rib_samples.phi_peri_sigma))
        pred_median = np.median(pred_phi_peri)
        pred_2sig =  np.percentile(pred_phi_peri, [2.5, 97.5])
        pred_sig = np.percentile(pred_phi_peri, [16, 84]) 
        ax[2].vlines(d['phi_peri'].values[i], pred_2sig[0], pred_2sig[1], lw=0.5, color=fmt['markeredgecolor'])
        ax[2].vlines(d['phi_peri'].values[i], pred_sig[0], pred_sig[1], lw=1, color=fmt['markeredgecolor'])
        ax[2].plot(d['phi_peri'].values[i], pred_median, **fmt) 

ax[2].plot([0, 0.15], [0, 0.15], '--', color=cor['primary_black'])
ax[3].plot([3, 10], [3, 10], '--', color=cor['primary_black'])
ax[0].set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax[0].set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
ax[1].set_xlabel('measured\nmembrane fraction', fontsize=6)
ax[1].set_ylabel('predicted\nmembrane fraction', fontsize=6)
ax[2].set_xlabel('measured\nperiplasmic fraction', fontsize=6)
ax[2].set_ylabel('predicted\nperiplasmic fraction', fontsize=6)
ax[3].set_xlabel('measured $S_A/V$ [µm$^{-1}$]', fontsize=6)
ax[3].set_ylabel('predicted $S_A/V$ [µm$^{-1}$]', fontsize=6)
plt.tight_layout()
plt.savefig('./plots/fig4_perturbations_predictions.pdf', bbox_inches='tight')