#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')
fit = pd.read_csv('../../data/mcmc/theory_growth_rate_prediction_summaries.csv')

fig, ax = plt.subplots(1, 1, figsize=(3,2))

inter_colors = {'95%':'pale_', '75%':'light_', '25%':'primary_', 'median':''}


for g, d in fit[fit['quantity']=='pred_lam_prot'].groupby('interval', sort=False):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=cor[f'{inter_colors[g]}red'], alpha=0.5)
for g, d in data.groupby('source'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['fg_protein_per_cell'], **fmt)

ax.legend()
ax.set_ylabel('$M_{prot}^{(tot)}$ [fg / cell]\ntotal protein', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

plt.savefig('../../figures/FigAX_protein_bayes.pdf')