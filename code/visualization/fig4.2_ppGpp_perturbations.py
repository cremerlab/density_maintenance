#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz
cor, pal = size.viz.matplotlib_style()
# Define constants
BETA = 1/0.4558

# Load the data, parameter samples, and wild-type predictions
data = pd.read_csv('../../data/compiled_measurements.csv')
par_samples = pd.read_csv('../../data/mcmc/fig3_inference_samples.csv')
quants = pd.read_csv('../../data/mcmc/fig3_fits.csv')

# Extract dependent quantities
phi_rib = data['phi_rib']
phi_mem = data['phi_mem']
phi_peri = data['phi_peri']

# Storage vectors for means and bounds
mean_vals = []
sig2_upper = [] 
sig2_lower = [] 
sig1_upper = [] 
sig1_lower = []

# Define the percentiles to compute
percs = [(2.5, 97.5), (16, 84)]

# Iterate through each measurement
for i in range(len(data)):
    # Compute the distribution of predicted SAV
    sav_pred = phi_mem[i] * par_samples['kappa'] / (2 * (1 + BETA * phi_rib[i] - phi_mem[i] - phi_peri[i]))

    # Compute the percentiles
    mean_vals.append(sav_pred.mean())
    for p, (low, high) in zip(percs, [(sig2_lower, sig2_upper), 
                                      (sig1_lower, sig1_upper)]):
        _p = np.percentile(sav_pred, p)     
        low.append(_p[0])
        high.append(_p[1])

# Update the dataframe
data['pred_sav_mean'] = mean_vals
data['1sig_lower'] = sig1_lower
data['1sig_upper'] = sig1_upper
data['2sig_lower'] = sig2_lower
data['2sig_upper'] = sig2_upper

#%% Define a graphical language for the perturbations
inducer_concs = {0:'pale_', 2:'primary_', 4:'', 100:''}
strain_colors = {'relA':  'red',
                 'meshI': 'gold'} 
markers = {'glucose': 'o', 'glucoseCAA': 's'}

#%% Plot the SAV theory and unceratainty with wildtype and perturbation data
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))

# Plot the predictions
sav_pred = quants[quants['quantity']=='theory_ppc']
ax.fill_between(sav_pred['phi_rib'], sav_pred['2sig_lower'], sav_pred['2sig_upper'],
                color=cor['pale_black'], alpha=0.5)
ax.fill_between(sav_pred['phi_rib'], sav_pred['1sig_lower'], sav_pred['1sig_upper'],
                color=cor['light_black'], alpha=0.5)
ax.plot(sav_pred['phi_rib'], sav_pred['mean_val'], '-', lw=1,
        color=cor['primary_black'], alpha=0.5)

# Plot the wildtype data
wt_data = data[data['strain']=='wildtype']
fmt = size.viz.style_point('This Study')
fmt['alpha'] = 0.5
ax.plot(wt_data['phi_rib'], wt_data['sav_inv_um'], **fmt)

# Plot the perturbation data
for g, d in data[data['strain']!='wildtype'].groupby(['strain', 'inducer_conc', 'carbon_source']):
    fmt = {'marker':markers[g[-1]],
           'markeredgecolor': cor[f'dark_{strain_colors[g[0]]}'],
           'markeredgewidth': 0.5,
           'markersize': 4,
           'markerfacecolor': cor[f'{inducer_concs[g[1]]}{strain_colors[g[0]]}'],
           'linestyle':'none'}

    ax.plot(d['phi_rib'], d['sav_inv_um'], **fmt)


# Add context
ax.set_xlim([0.1, 0.25])
ax.set_ylim([4, 8])
ax.set_xlabel('ribosomal proteome fraction\n$\phi_{rib}$', fontsize=6)
ax.set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)

#%% 
fig, ax = plt.subplots(1, 2, figsize=(2.5, 1.5))