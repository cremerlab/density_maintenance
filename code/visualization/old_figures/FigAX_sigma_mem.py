#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
fit = pd.read_csv('../../data/mcmc/theory_growth_rate_prediction_summaries.csv')

fig, ax = plt.subplots(1, 1, figsize=(3,2))
ax.set_ylim([0, 6])
inter_colors = {'95%':'pale_', '75%':'light_', '25%':'primary_', 'median':''}

for g, d in fit[fit['quantity']=='rho_mem_pred'].groupby('interval', sort=False):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=cor[f'{inter_colors[g]}red'], alpha=0.5)
for g, d in data[data['quantity']=='ms_rho_mem'].groupby('source'):
    fmt = size.viz.style_point(g, alpha=0.5)
    ax.plot(d['growth_rate_hr'], d['median_value'], **fmt,ms=5)

ax.legend(fontsize=5, bbox_to_anchor=(1,1))
ax.set_ylabel('$\sigma_{mem}$ [fg / µm$^2$\naverage membrane protein areal density', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

plt.savefig('../../figures/FigAX_sigma_mem_bayes.pdf')