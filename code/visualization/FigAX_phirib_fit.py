#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
data = data[(data['source']!='Si et al. 2017') & (data['source']!='Taher-Araghi et al. 2015')]
fit = pd.read_csv('../../data/mcmc/theory_growth_rate_prediction_summaries.csv')

fig, ax = plt.subplots(1, 1, figsize=(3,2))

inter_colors = {'95%':'pale_', '75%':'light_', '25%':'primary_', 'median':''}

for g, d in fit[fit['quantity']=='phi_Rb_pred'].groupby('interval', sort=False):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=cor[f'{inter_colors[g]}red'], alpha=0.5, zorder=1000)
for g, d in data.groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    ax.plot(d['growth_rate_hr'], d['mass_fraction'], **fmt,ms=5)

ax.legend(fontsize=5, bbox_to_anchor=(1,1))
ax.set_ylabel('$\phi_{rib}$\nribosomal proteome fraction', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

plt.savefig('../../figures/FigAX_phi_rib_bayes.pdf')