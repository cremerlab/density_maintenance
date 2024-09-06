#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/literature/collated_mass_fractions_empirics.csv')
data = data[data['localization']=='periplasm']
fit = pd.read_csv('../../data/mcmc/theory_growth_rate_prediction_summaries.csv')

fig, ax = plt.subplots(1, 1, figsize=(3,2))

inter_colors = {'95%':'pale_', '75%':'light_', '25%':'primary_', 'median':''}

for g, d in fit[fit['quantity']=='phi_peri_lam_pred'].groupby('interval', sort=False):
    ax.fill_between(d['growth_rate_hr'], d['lower'], d['upper'], color=cor[f'{inter_colors[g]}red'], alpha=0.5)

for g, d in data.groupby('dataset_name'):
    fmt = size.viz.style_point(g)
    ax.plot(d['growth_rate_hr'], d['mass_frac'], **fmt)
ax.set_ylim([0, 0.12])
ax.legend(fontsize=5, loc='upper right')
ax.set_ylabel('$\phi_{peri}$ \nperiplasmic proteome fraction', fontsize=6)
ax.set_xlabel('growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

plt.savefig('../../figures/FigAX_phi_peri_bayes.pdf')