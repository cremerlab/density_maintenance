#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

size_emp = pd.read_csv('../../data/mcmc/size_data_empirical_summaries_wide.csv')
wt_data = pd.read_csv('../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
preds = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')


fig, ax = plt.subplots(2, 2, figsize=(4,4), sharex=True)

for i in range(2):
    ax[i,0].set_ylabel('$S_A/V$ [µm$^{-1}$]\nsurface-to-volume', fontsize=6)
    ax[i,1].set_ylabel('$w$ [µm]\naverage width', fontsize=6)
    ax[i, 0].set_ylim(3.5, 9)
    ax[i, 1].set_ylim(0.45, 1.2)
    ax[1, i].set_xlabel('RNA-to-protein\n$M_{RNA}/M_{prot}^{(tot)}$', fontsize=6)
for a in ax.ravel():
    a.set_xlim([0.1, 0.7])

for g, d in size_emp.groupby('source'):
    fmt = size.viz.style_point(g)
    for i in range(2): 
        for j, _d in enumerate([d['surface_to_volume'], d['width']]):
            ax[i, j].hlines(_d, d['97.5%']/0.4558, d['2.5%']/0.4558, lw=0.5, color=fmt['color'])
            ax[i, j].plot(d['median_value'] / 0.4558, _d, **fmt)
            
for g, d in wt_data.groupby('carbon_source'):
    phi = d[d['quantity']=='phiRb']
    width = d[d['quantity']=='width']
    sav = d[d['quantity']=='surface_to_volume']
    for j, _d in enumerate([sav, width]):
        for i in range(2):
            ax[i, j].hlines(_d['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558, lw=1, color=cor['primary_black'])
            ax[i, j].vlines(phi['median_value']/0.4558, _d['2.5%'], _d['97.5%'], lw=1, color=cor['primary_black'])
            ax[i, j].plot(phi['median_value']/0.4558, _d['median_value'], 'o', markerfacecolor='w', markeredgecolor=cor['primary_black'], markeredgewidth=1)


inter_colors = {'95%':'pale_', '75%':'light_', '25%':'primary_', 'median':''}
c = ['blue', 'blue', 'green', 'green']
ax = ax.ravel()
for i, q in enumerate(['SAV_theory_simple', 'width_theory_simple', 'SAV_theory', 'width_theory']):
    d = preds[preds['quantity']==q]
    for g, _d in d.groupby('interval', sort=False):
        ax[i].fill_between(_d['phiRb']/0.4558, _d['lower'], _d['upper'], color=cor[f'{inter_colors[g]}{c[i]}'], alpha=0.75)
plt.tight_layout()