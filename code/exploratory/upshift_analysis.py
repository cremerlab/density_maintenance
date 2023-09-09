# %%
import numpy as np
import pandas as pd
import size.viz
import scipy.stats
import size.fluxparity
import matplotlib.pyplot as plt
const = size.fluxparity.load_constants()
cor, pal = size.viz.matplotlib_style()

# %%
data = pd.read_csv('../../data/literature/woldringh1980_upshift.csv')
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')

lam_pre = np.log(2) / (72/60)  # From caption of figure 1
lam_post = np.log(2) / (24/60)  # From caption of figure 2

# Estimate the steady-state ribosomal allocation for these growth rates
popt = scipy.stats.linregress(
    phiRb_data['growth_rate_hr'], phiRb_data['mass_fraction'])
phiRb_pre = popt[1] + popt[0] * lam_pre
phiRb_pre = 0.05
phiRb_post = popt[1] + popt[0] * lam_post

nu_pre = size.fluxparity.estimate_nu_FPM(
    phiRb_pre, lam_pre, const, nu_buffer=10, phi_O=0.55 + 0.12 + 0.08)
nu_post = size.fluxparity.estimate_nu_FPM(
    phiRb_post, lam_post, const, nu_buffer=10, phi_O=0.55 + 0.12 + 0.02)

# %%
_args = [{'gamma_max': const['gamma_max'],
          'nu_max':nu_pre,
          'Kd_TAA': const['Kd_TAA'],
          'Kd_TAA_star': const['Kd_TAA_star'],
          'kappa_max': const['kappa_max'],
          'phi_O': 0.65,
          'phi_peri': 0.1,
          'phi_mem': 0.12,
          'tau':const['tau']},

         {'gamma_max': const['gamma_max'],
          'nu_max':nu_post,
          'Kd_TAA': const['Kd_TAA'],
          'Kd_TAA_star': const['Kd_TAA_star'],
          'kappa_max': const['kappa_max'],
          'phi_O': 0.55,
          'phi_peri': 0.12,
          'phi_peri': 0.02,
          'tau':const['tau']}
         ]
shift_df = size.fluxparity.nutrient_shift_FPM(_args, 1, 6)
# %%
fig, ax = plt.subplots(1, 2, figsize=(4, 2))
ax[0].set_ylabel('ribosomal mass fraction', fontsize=6)
ax[0].set_xlabel('time [hr]', fontsize=6)
ax[0].plot(shift_df['shifted_time'].values, shift_df['M_Rb'] / shift_df['M'], '-',
           color=cor['primary_blue'], lw=1.5)
ax[1].plot(data['time_hr'], data['width_um'], 'o')
phi_peri = 0.02
phi_mem = 0.12
kappa = 109
alpha = 3.58
preshift = shift_df[shift_df['shifted_time'] < 0]
postshift = shift_df[shift_df['shifted_time'] >= 0]
theo_pre = (24 * alpha / (3 * alpha - 1)) * (1 +
                                             preshift['M_Rb']/(0.4558 * preshift['M']) - 0.1 - phi_mem) / (kappa * phi_mem)
theo_post = (24 * alpha / (3 * alpha - 1)) * (1 + postshift['M_Rb']/(
    0.4558 * postshift['M']) - 0.02 - phi_mem) / (kappa * phi_mem)
ax[1].plot(preshift['shifted_time'], theo_pre, 'k-')
ax[1].plot(postshift['shifted_time'], theo_post, 'k-')
