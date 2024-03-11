#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
wt_post = pd.read_csv('../../data/mcmc/wildtype_posterior_parameter_summaries.csv')
temp_post = pd.read_csv('../../data/mcmc/temperature_perturbation_parameter_summaries.csv')
theo_post = pd.read_csv('../../data/mcmc/theory_phiRb_prediction_summaries.csv')

c_mapper = {'glycerol': 's', 'glucoseCAA': 'v'}
t_mapper = {25:cor['primary_blue'], 30:cor['light_blue'],
            42:cor['primary_red']}



creds = {'95%': cor['pale_green'], '75%': cor['light_green'], '25%': cor['primary_green'], 'median':cor['green']}

fig, ax = plt.subplots(2, 2, figsize=(3, 3))
ax[0,0].set_xlabel(r'1/T [$\times$10$^{3}$ C$^{-1}$]', fontsize=6)
ax[0,0].set_ylabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0,0].set_xlim([20, 45])
ax[0,0].set_ylim([0.15, 1.5])
ax[0, 1].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0, 1].set_ylabel('$M_{RNA}^{(tot)}/M_{prot}^{(tot)}$' + '\n RNA-to-protein', fontsize=6)
ax[1,0].set_xlim([0.08, 0.5])
ax[1,1].set_xlim([0.08, 0.5])
ax[1,0].set_ylim([5, 9])
ax[1, 1].set_ylim([0.45, 0.9])
for i in range(2):
    ax[1, i].set_xlabel('$M_{RNA}^{(tot)}/M_{prot}^{(tot)}$' + '\n RNA-to-protein', fontsize=6)
ax[1,0].set_ylabel('$S_A/V$\nsurface-to-volume [µm$^{-1}$]', fontsize=6)
ax[1,1].set_ylabel('w\naverage width [µm]', fontsize=6)

_ax = {'SAV_theory':ax[1, 0], 'width_theory':ax[1,1]}
for g, d in theo_post[theo_post['quantity'].isin(['SAV_theory', 'width_theory'])].groupby(['quantity','interval'], sort=False):
    _ax[g[0]].fill_between(d['phiRb']/0.4558, d['lower'], d['upper'], color=creds[g[1]], alpha=0.5)

for c, s in c_mapper.items():
    lam = wt_post[(wt_post['carbon_source']==c) & 
                  (wt_post['overexpression']=='none') & 
                  (wt_post['quantity'] == 'growth_rate')]
    
    phi = wt_post[(wt_post['carbon_source']==c) & 
                  (wt_post['overexpression']=='none') & 
                  (wt_post['quantity'] == 'phiRb')]

    sav = wt_post[(wt_post['carbon_source']==c) & 
                  (wt_post['overexpression']=='none') & 
                  (wt_post['quantity'] == 'surface_to_volume')]

    width = wt_post[(wt_post['carbon_source']==c) & 
                  (wt_post['overexpression']=='none') & 
                  (wt_post['quantity'] == 'width')]
 
    ax[0,0].vlines(27.01, lam['2.5%'], lam['97.5%'], lw=1, color=cor['primary_black'])
    ax[0,0].plot(27.02, lam['median_value'], s, markerfacecolor='w', 
                markeredgecolor=cor['primary_black'], markeredgewidth=1, ms=4)
    ax[0,1].vlines(lam['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558,
                 lw=1, color=cor['primary_black'])
    ax[0,1].hlines(phi['median_value']/0.4558, lam['2.5%'], lam['97.5%'],
                 lw=1, color=cor['primary_black'])
    ax[0,1].plot(lam['median_value'], phi['median_value']/0.4558, s, markerfacecolor='w',
            markeredgecolor=cor['primary_black'], markeredgewidth=1, ms=4)

    ax[1,0].vlines(phi['median_value']/0.4558, sav['2.5%'], sav['97.5%'],
                 lw=1, color=cor['primary_black'])
    ax[1,0].hlines(sav['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558,
                 lw=1, color=cor['primary_black'])
    ax[1,0].plot(phi['median_value']/0.4558, sav['median_value'], s, markerfacecolor='w',
            markeredgecolor=cor['primary_black'], markeredgewidth=1, ms=4)

    ax[1,1].vlines(phi['median_value']/0.4558, width['2.5%'], width['97.5%'],
                 lw=1, color=cor['primary_black'])
    ax[1,1].hlines(width['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558,
                 lw=1, color=cor['primary_black'])
    ax[1,1].plot(phi['median_value']/0.4558, width['median_value'], s, markerfacecolor='w',
            markeredgecolor=cor['primary_black'], markeredgewidth=1, ms=4)

for g, d in temp_post.groupby(['carbon_source', 'temperature']):
    lam = d[d['quantity']=='growth_rate']
    width = d[d['quantity']=='width']
    sav = d[d['quantity']=='surface_to_volume']
    phi = d[d['quantity']=='phiRb']
    x = 1000 / g[1]

    ax[0,0].vlines(x, lam['2.5%'], lam['97.5%'], lw=1, color=t_mapper[g[1]])
    ax[0,0].plot(x, lam['median_value'], c_mapper[g[0]], markerfacecolor='w',
            markeredgecolor=t_mapper[g[1]], markeredgewidth=1, ms=4)
    ax[0,1].vlines(lam['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558,
                 lw=1, color=t_mapper[g[1]])
    ax[0,1].hlines(phi['median_value']/0.4558, lam['2.5%'], lam['97.5%'],
                 lw=1, color=t_mapper[g[1]])
    ax[0,1].plot(lam['median_value'], phi['median_value']/0.4558, c_mapper[g[0]], markerfacecolor='w',
            markeredgecolor=t_mapper[g[1]], markeredgewidth=1, ms=4)

    ax[1,0].vlines(phi['median_value']/0.4558, sav['2.5%'], sav['97.5%'],
                 lw=1, color=t_mapper[g[1]])
    ax[1,0].hlines(sav['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558,
                 lw=1, color=t_mapper[g[1]])
    ax[1,0].plot(phi['median_value']/0.4558, sav['median_value'], c_mapper[g[0]], markerfacecolor='w',
            markeredgecolor=t_mapper[g[1]], markeredgewidth=1, ms=4)

    ax[1,1].vlines(phi['median_value']/0.4558, width['2.5%'], width['97.5%'],
                 lw=1, color=t_mapper[g[1]])
    ax[1,1].hlines(width['median_value'], phi['2.5%']/0.4558, phi['97.5%']/0.4558,
                 lw=1, color=t_mapper[g[1]])
    ax[1,1].plot(phi['median_value']/0.4558, width['median_value'], c_mapper[g[0]], markerfacecolor='w',
            markeredgecolor=t_mapper[g[1]], markeredgewidth=1, ms=4)

plt.subplots_adjust(hspace=0.3, wspace=0.5)
plt.savefig('../../figures/FigSX_temperature_growth_rate.pdf', bbox_inches='tight')