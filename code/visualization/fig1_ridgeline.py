# %%
import size.viz
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import KernelDensity
cor, pal = size.viz.matplotlib_style()
gr_data = pd.read_csv('../../data/mcmc/wildtype_growth_rate_summary.csv')
size_data = pd.read_csv(
    '../processing/microscopy/wildtype_size_measurement/output/wildtype_size_measurements.csv')
size_data = size_data[size_data['carbon_source'] != 'ezMOPS']
size_data = size_data[(size_data['width_median'] >= 0.2)
                      & (size_data['length'] <= 6.5)]
prot_data = pd.read_csv(
    '../../data/protein_quantification/wildtype_simple_quantification.csv')
cmap = sns.color_palette("ch:start=.2,rot=-.3", n_colors=7).as_hex()
carbons = ['LB', 'glucoseCAA', 'glucose', 'glycerol', 'sorbitol', 'acetate']
carb_cor = {c: cmap[-i - 1] for i, c in enumerate(carbons)}

# %%
_carbons = ['LB', 'glucose +\nCAA', 'glucose',
            'glycerol', 'sorbitol', 'acetate']
# Select the carbon sources
carbons = ['LB', 'glucoseCAA', 'glucose', 'glycerol', 'sorbitol', 'acetate']
data = size_data[size_data['carbon_source'].isin(carbons)]
ind = {k: i for i, k in enumerate(carbons)}
wbins = np.linspace(0.1, 2, 50)
lbins = np.linspace(0.5, 7, 25)

nudge = 0.1

fig, ax = plt.subplots(1, 2, figsize=(3, 1.5), sharey=True)
ax[0].set_xlabel('cell width [µm]', fontsize=6)
ax[1].set_xlabel('cell length [µm]', fontsize=6)
for a in ax:
    a.grid(False)
    a.set_facecolor('white')
    a.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
ax[0].set_yticklabels(_carbons, fontweight='bold', fontsize=6)

for g, d in size_data.groupby(['carbon_source']):
    whist, _ = np.histogram(d['width_median'], bins=wbins, density=True)
    lhist, _ = np.histogram(d['length'], bins=lbins, density=True)
    idx = 0
    for bins, hist in zip([wbins, lbins], [whist, lhist]):
        hist *= hist.max()**-1
        hist *= 0.2
        ax[idx].step(bins[:-1], hist + nudge * ind[g], color=carb_cor[g])
        ax[idx].fill_between(bins[:-1], nudge * ind[g], hist + nudge * ind[g],
                             alpha=0.5, step='pre',
                             zorder=len(ind) - ind[g], color=carb_cor[g])
        idx += 1
ax[1].set_xlim([0.5, 7])
plt.tight_layout()
plt.savefig('../../figures/fig1_data_size_histograms.pdf', bbox_inches='tight')


# %%
#  Plot the means
fig, ax = plt.subplots(1, 2, figsize=(3, 1), sharex=True)
ax[0].set_ylim([0.6, 1.1])
ax[1].set_ylim([1.75, 4])
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('average\nwidth [µm]', fontsize=6)
ax[1].set_ylabel('average\nlength [µm]', fontsize=6)
ax[0].set_yticks([0.6, 0.8, 1.0, 1.2])
ax[1].set_yticks([2, 3, 4])
for g, d in size_data.groupby(['carbon_source']):
    gr = gr_data[gr_data['carbon_source'] == g]['median'].values[0]
    means_grouped = d.groupby(['date']).agg('mean').reset_index()
    mean_w = means_grouped['width_median'].mean()
    mean_l = means_grouped['length'].mean()
    sem_w = means_grouped['width_median'].std() / np.sqrt(len(means_grouped))
    sem_l = means_grouped['length'].std() / np.sqrt(len(means_grouped))
    ax[0].errorbar(gr, mean_w, yerr=sem_w, marker='o',
                   linewidth=1,  ms=4, color=carb_cor[g])
    ax[1].errorbar(gr, mean_l, yerr=sem_l, marker='o',
                   linewidth=1,  ms=4, color=carb_cor[g])
plt.tight_layout()
plt.savefig('../../figures/fig1_avg_length_width_gr.pdf', bbox_inches='tight')

# %%
fig, ax = plt.subplots(1, 2, figsize=(3.5, 2), sharex=True)
for a in ax:
    a.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('surface area to  volume ratio [µm$^{-1}$]', fontsize=6)
ax[1].set_ylabel(
    'periplasmic protein per biomass\n[µg$\cdot$mL$^{-1}\cdot$OD$_{600nm}^{-1}$]', fontsize=6)
ax[0].set_ylim([5, 7])
ax[1].set_ylim([5, 30])
ax[0].set_yticks([5, 5.5, 6, 6.5, 7])


# Plot the surface area to volume
for g, d in size_data.groupby(['carbon_source']):
    if g == 'LB':
        continue
    gr = gr_data[gr_data['carbon_source'] == g]['median'].values[0]
    means_grouped = d.groupby(['date']).agg('mean').reset_index()
    mean_sav = means_grouped['surface_to_volume'].mean()
    sem_sav = means_grouped['surface_to_volume'].std() / \
        np.sqrt(len(means_grouped))
    ax[0].errorbar(gr, mean_sav, sem_sav, marker='o',
                   ms=4, color=carb_cor[g], lw=1)

for g, d in prot_data.groupby(['carbon_source']):
    if g == 'LB':
        continue
    gr = gr_data[gr_data['carbon_source'] == g]['median'].values[0]
    ax[1].errorbar(gr, d['m_peri_per_biomass'].mean(), d['m_peri_per_biomass'].std() / np.sqrt(len(d)), color=carb_cor[g],
                   marker='o', ms=4, lw=1)
plt.tight_layout()
plt.savefig('../../figures/fig1_SVR_Mperi.pdf', bbox_inches='tight')
