# %%
import size.viz
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import KernelDensity
cor, pal = size.viz.matplotlib_style()
size_data = pd.read_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
size_data = size_data[size_data['carbon_source'] != 'ezMOPS']
gr_data = pd.read_csv('../../data/mcmc/wildtype_growth_rate_summary.csv')
cmap = sns.color_palette("ch:start=.2,rot=-.3", n_colors=7).as_hex()
carbons = ['LB', 'glucoseCAA', 'glucose', 'glycerol', 'sorbitol', 'acetate']
carb_cor = {c: cmap[-i - 1] for i, c in enumerate(carbons)}

# %%
fig, ax = plt.subplots(2, 1, figsize=(2.5, 2), sharex=True)
ax[0].set_ylabel('cell length [µm]', fontsize=6)
ax[1].set_ylabel('cell width [µm]', fontsize=6)
ax[1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)

for g, d in size_data.groupby(['carbon_source']):
    _len = d[d['parameter'] == 'length_um']
    _width = d[d['parameter'] == 'width_um']
    _gr = gr_data[(gr_data['carbon_source'] == g)]
    if (len(_gr) == 0) | (len(_len) == 0):
        continue
    if g == 'glucoseCAA':
        continue
    # Plot the error bars for length
    ax[0].vlines(_gr['median'], _len["97.5%"], _len["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[0].vlines(_gr['median'], _len["87.5%"], _len["12.5%"],  lw=1,
                 color=carb_cor[g])
    ax[0].hlines(_len['median'], _gr["97.5%"], _gr["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[0].hlines(_len['median'], _gr["87.5%"], _gr["12.5%"], lw=1,
                 color=carb_cor[g])
    ax[0].plot(_gr['median'], _len['median'], 'o', markerfacecolor='w',
               markeredgecolor=carb_cor[g], ms=2, markeredgewidth=0.75)

    # Plot the error bars for width
    ax[1].vlines(_gr['median'], _width["97.5%"], _width["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[1].vlines(_gr['median'], _width["87.5%"], _width["12.5%"],  lw=1,
                 color=carb_cor[g])
    ax[1].hlines(_width['median'], _gr["97.5%"], _gr["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[1].hlines(_width['median'], _gr["87.5%"], _gr["12.5%"], lw=1,
                 color=carb_cor[g])
    ax[1].plot(_gr['median'], _width['median'], 'o', markerfacecolor='w',
               markeredgecolor=carb_cor[g], ms=2, markeredgewidth=0.75)
plt.tight_layout()
plt.savefig('../../Fig1_cell_length_width.pdf')

# %%
for g, d in size_data.groupby(['carbon_source']):
    _len = d[d['parameter'] == 'length_um']
    _width = d[d['parameter'] == 'width_um']
    _gr = gr_data[(gr_data['carbon_source'] == g)]
    if g == 'glucoseCAA':
        continue
    # Plot the error bars for length
    ax[0].vlines(_gr['median'], _len["97.5%"], _len["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[0].vlines(_gr['median'], _len["87.5%"], _len["12.5%"],  lw=1,
                 color=carb_cor[g])
    ax[0].hlines(_len['median'], _gr["97.5%"], _gr["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[0].hlines(_len['median'], _gr["87.5%"], _gr["12.5%"], lw=1,
                 color=carb_cor[g])
    ax[0].plot(_gr['median'], _len['median'], 'o', markerfacecolor='w',
               markeredgecolor=carb_cor[g], ms=2, markeredgewidth=0.75)

    # Plot the error bars for width
    ax[1].vlines(_gr['median'], _width["97.5%"], _width["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[1].vlines(_gr['median'], _width["87.5%"], _width["12.5%"],  lw=1,
                 color=carb_cor[g])
    ax[1].hlines(_width['median'], _gr["97.5%"], _gr["2.5%"], lw=0.5,
                 color=carb_cor[g])
    ax[1].hlines(_width['median'], _gr["87.5%"], _gr["12.5%"], lw=1,
                 color=carb_cor[g])
    ax[1].plot(_gr['median'], _width['median'], 'o', markerfacecolor='w',
               markeredgecolor=carb_cor[g], ms=2, markeredgewidth=0.75)
plt.tight_layout()
plt.savefig('../../Fig1_cell_length_width.pdf')


# %%
# Load the sampling
data = pd.read_csv('../../data/mcmc/wildtype_hyperparameter_size_samples.csv')
data = pd.read_csv('../../data/')

# Select the carbon sources
carbons = ['LB', 'glucoseCAA', 'glucose', 'glycerol', 'sorbitol', 'acetate']
_carbons = ['LB', 'glucose +\nCAA', 'glucose',
            'glycerol', 'sorbitol', 'acetate']
data = data[data['carbon_source'].isin(carbons)]


# Set up the figure axes
fig, ax = plt.subplots(1, 2, figsize=(2.8, 1.8), sharey=True)
ax[0].set_xlabel('cell length [µm]', fontsize=6)
ax[0].set_xticks([1, 2, 3, 4])
ax[1].set_xlabel('cell width [µm]', fontsize=6)
ind = {k: i for i, k in enumerate(carbons)}


w_range = np.linspace(0.1, 1.5, 200)[:, np.newaxis]
l_range = np.linspace(1.5, 4, 200)[:, np.newaxis]
nudge = 0.6
for g, d in data.groupby(['carbon_source']):
    print(g)
    if g == 'glucoseCAA':
        continue
    width_kde = KernelDensity(kernel='gaussian', bandwidth=0.1).fit(
        d['width_um'].values[:, np.newaxis])
    width_dens = np.exp(width_kde.score_samples(w_range))
    width_dens *= width_dens.max()**-1
    len_kde = KernelDensity(kernel='gaussian', bandwidth=0.1).fit(
        d['length_um'].values[:, np.newaxis])
    len_dens = np.exp(len_kde.score_samples(l_range))
    len_dens *= len_dens.max()**-1
    ax[1].plot(w_range[:, 0], width_dens + nudge *
               ind[g], zorder=-ind[g], color=carb_cor[g])
    ax[1].fill_between(w_range[:, 0],  nudge * ind[g] + width_dens.min(),
                       width_dens + nudge * ind[g], alpha=0.5, zorder=-ind[g], color=carb_cor[g])
    ax[0].plot(l_range[:, 0], len_dens + nudge * ind[g], color=carb_cor[g])
    ax[1].fill_between(w_range[:, 0],  nudge * ind[g] + width_dens.min(),
                       width_dens + nudge * ind[g], alpha=0.5, zorder=-ind[g], color=carb_cor[g])
    ax[0].fill_between(l_range[:, 0],  nudge * ind[g] + len_dens.min(),
                       len_dens + nudge * ind[g], alpha=0.5, zorder=-ind[g], color=carb_cor[g])

for i, a in enumerate(ax):
    if i == 0:
        a.set_yticks(nudge * np.arange(len(carbons)))
        a.set_yticklabels(_carbons)
    a.set_facecolor('none')
    a.grid(False)
plt.tight_layout()
plt.savefig('../../figures/fig1_cell_length_width_posteriors.pdf')
# %%
_kde = KernelDensity(kernel='gaussian', bandwidth=0.1).fit(
    d['width_um'].values[:, np.newaxis])
log_dens = _kde.score_samples(w_range)
# plt.hist(d['width_um'])
plt.plot(w_range[:, 0], np.exp(log_dens))

# %%
fig, ax = plt.subplots(1, 2, figsize=(2.8, 1.8), sharey=True)
ax[0].set_xlabel('cell length [µm]', fontsize=6)
ax[0].set_xticks([1, 2, 3, 4])
ax[1].set_xlabel('cell width [µm]', fontsize=6)
ind = {k: i for i, k in enumerate(carbons)}


w_range = np.linspace(0.1, 1.5, 200)[:, np.newaxis]
l_range = np.linspace(1.5, 4, 200)[:, np.newaxis]
nudge = 0.6
width_bins = np.linspace(0.1, 2, 400)
length_bins = np.linspace(1, 4, 400)
for g, d in data.groupby(['carbon_source']):
    print(g)
    # if g == 'glucoseCAA':
    # continue
    bins, hist = np.histogram(d['width_um'])
    ax[1].plot(w_range[:, 0], width_dens + nudge *
               ind[g], zorder=-ind[g], color=carb_cor[g])
    ax[1].fill_between(w_range[:, 0],  nudge * ind[g] + width_dens.min(),
                       width_dens + nudge * ind[g], alpha=0.5, zorder=-ind[g], color=carb_cor[g])
    ax[0].plot(l_range[:, 0], len_dens + nudge * ind[g], color=carb_cor[g])
    ax[1].fill_between(w_range[:, 0],  nudge * ind[g] + width_dens.min(),
                       width_dens + nudge * ind[g], alpha=0.5, zorder=-ind[g], color=carb_cor[g])
    ax[0].fill_between(l_range[:, 0],  nudge * ind[g] + len_dens.min(),
                       len_dens + nudge * ind[g], alpha=0.5, zorder=-ind[g], color=carb_cor[g])

for i, a in enumerate(ax):
    if i == 0:
        a.set_yticks(nudge * np.arange(len(carbons)))
        a.set_yticklabels(_carbons)
    a.set_facecolor('none')
    a.grid(False)
plt.tight_layout()
plt.savefig('../../figures/fig1_cell_length_width_posteriors.pdf')
