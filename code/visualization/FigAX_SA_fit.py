#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/literature/collated_literature_size_data.csv')
fit = pd.read_csv('../../data/mcmc/theory_growth_rate_prediction_summaries.csv')

fig, ax = plt.subplots(1, 1, figsize=(3,2))

inter_colors = {'95%':'pale_', '75%':'light_', '25%':'primary_', 'median':''}

for g, d in fit[fit['quantity']=='SA_fit'].groupby('interval', sort=False):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=cor[f'{inter_colors[g]}red'], alpha=0.5)

for g, d in data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['surface_area_um2'], **fmt)
ax.set_ylim([0, 15])
ax.legend(fontsize=6, loc='upper left')
ax.set_ylabel('$S_A$ [µm$^2$]\naverage cell surface area', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

plt.savefig('../../figures/FigAX_sa_bayes.pdf')