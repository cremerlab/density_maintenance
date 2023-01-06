# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import seaborn as sns
cor, pal = size.viz.matplotlib_style()
cmap = sns.color_palette("ch:start=.2,rot=-.3", n_colors=7).as_hex()
carb_cor = {k: v for k, v in zip(
    ['acetate', 'sorbitol', 'glycerol', 'glucose', 'glucoseCAA'], cmap)}

voldata = pd.read_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
protdata = pd.read_csv(
    '../../data/protein_quantification/mcmc/wildtype_protein_per_biomass_hyperparameter_summary.csv')
delta = 0.025

theta_range = np.linspace(5, 30)
theory = theta_range * (100 * delta)**-1

fig, ax = plt.subplots(1, 1, figsize=(2, 2))

for g, d in protdata.groupby(['carbon_source']):
    _sav = voldata[(voldata['carbon_source'] == g) &
                   (voldata['parameter'] == 'SAV_inv_um')]
    # ax.vlines(d['median'], _sav['2.5%'],
    #   _sav['97.5%'], lw=0.5, color=carb_cor[g])
    # ax.vlines(d['median'], _sav['12.5%'],
    #   _sav['87.5%'], lw=1, color=carb_cor[g])
    # ax.hlines(_sav['median'], d['2.5%'], d['97.5%'], lw=0.5, color=carb_cor[g])
    # ax.hlines(_sav['median'], d['12.5%'], d['87.5%'], lw=1, color=carb_cor[g])
    ax.plot(d['median']/2, _sav['median'], 'o', color=carb_cor[g], ms=4)

# ax.set_xlim([20, 50])
# ax.set_ylim([4, 7])
ax.set_xlabel('periplasmic protein per biomass')
ax.set_ylabel('SAV')
ax.plot(theta_range, theory, 'k--')
