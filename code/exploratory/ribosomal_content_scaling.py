# %%
import numpy as np
import pandas as pd
import size.viz
import growth.model
import matplotlib.pyplot as plt
cor, pal = size.viz.matplotlib_style()
mapper = size.viz.lit_mapper()
const = growth.model.load_constants()

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
rib = ms_data[ms_data['localization'] == 'ribosomal sector']


nu_range = np.linspace(0.1, 20, 200)
opt_phiRb = growth.model.phiRb_optimal_allocation(
    const['gamma_max'], nu_range, const['Kd_cpc'], const['phi_O'])
lam = growth.model.steady_state_growth_rate(
    const['gamma_max'], opt_phiRb, nu_range, const['Kd_cpc'], const['phi_O'])

# model pred
k = 150
beta = 0.4558
alpha = 4
phi_mem = 0.1
pred = (24 * alpha / (k * (3 * alpha - 1))) * (1 + beta * opt_phiRb)/phi_mem

fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))
for g, d in rib.groupby(['dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['mass_frac'], mapper[g]['m'],
               markeredgewidth=0.25, markeredgecolor=cor['primary_black'],
               color=mapper[g]['c'], alpha=0.5, label=g)
    ax[1].plot(d['mass_frac'], d['width'], mapper[g]['m'],
               markeredgewidth=0.25, markeredgecolor=cor['primary_black'],
               color=mapper[g]['c'], alpha=0.5, label=g)

ax[0].plot(lam, opt_phiRb, '-', color=cor['primary_blue'],
           lw=2, label='optimal allocation theory')
# ax[1].plot(opt_phiRb, '-', color=cor['primary_blue'], lw=2)

ax[1].plot(opt_phiRb, pred, '-', color=cor['primary_blue'],
           lw=2, label='density maintenance theory')
ax[0].set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax[0].set_ylabel('ribosomal mass fraction $\phi_{Rb}$', fontsize=6)
ax[1].set_xlabel('ribosomal mass fraction $\phi_{Rb}$', fontsize=6)
ax[1].set_ylabel('average cell width [µm]', fontsize=6)
ax[0].legend(fontsize=5)
ax[1].legend(fontsize=5)
