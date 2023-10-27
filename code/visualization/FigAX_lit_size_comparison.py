#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

size_data = pd.read_csv('../../data/literature/collated_literature_size_data.csv')
size_data = size_data[~size_data['source'].isin(['Si et al. 2017', 'Taher-Araghi et al. 2015', 'Basan et al. 2015'])]
wt_data = pd.read_csv('../../data/mcmc/wildtype_posterior_parameter_summaries.csv')

fig, ax = plt.subplots(1, 2, figsize=(4, 1.5))
ax[0].set_xlim([0, 2.2])
ax[0].set_ylim([0.45, 1.2])
ax[1].set_ylim([1, 5])
for a in ax:
    a.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('$w$ [µm]\naverage cell width', fontsize=6)
ax[1].set_ylabel('$\ell$ [µm]\naverage cell length', fontsize=6)

for g, d in size_data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax[0].plot(d['growth_rate_hr'], d['width_um'], **fmt)
    ax[1].plot(d['growth_rate_hr'], d['length_um'], **fmt)


x = 0
for g, d in wt_data[wt_data['quantity'].isin(['width', 'length', 'growth_rate'])].groupby('carbon_source'):
    lam = d[d['quantity']=='growth_rate'] 
    w = d[d['quantity']=='width']
    l = d[d['quantity']=='length']
    if x == 0:
        label = 'This study'
    else:
        label = '__nolegend__'
    for i, p in enumerate([w, l]):
        ax[i].vlines(lam['median_value'], p['2.5%'], p['97.5%'], lw=1, color=cor['primary_black'], label='__nolegend__')
        ax[i].hlines(p['median_value'], lam['2.5%'], lam['97.5%'], lw=1, color=cor['primary_black'], label='__nolegend__')
        ax[i].plot(lam['median_value'], p['median_value'], 'o', ms=5, markeredgecolor=cor['primary_black'],
                   markeredgewidth=1, markerfacecolor='w', label=label)
        x += 1
ax[1].legend(fontsize=5, bbox_to_anchor=(1,1))
plt.savefig('../../figures/FigAX_lit_size_comparison.pdf')