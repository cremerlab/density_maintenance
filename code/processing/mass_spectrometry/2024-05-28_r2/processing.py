#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
import scipy.stats
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('./raw/2024-05-28_biuret_calibration.csv')

popt = scipy.stats.linregress(data['nominal_BSA_conc_ug_mL'], data['od_555nm'])

conc_range = np.linspace(0, 1100, 100)
fit = popt[1] + popt[0] * conc_range

fig, ax = plt.subplots(1,1, figsize=(2,1))
ax.plot(data['nominal_BSA_conc_ug_mL'], data['od_555nm'], 'o', color=cor['primary_black'],
        ms=5)
ax.plot(conc_range, fit, '-', color=cor['primary_red'], lw=1)
ax.set_xlabel('measured BSA concentration [ug/mL]', fontsize=6)
ax.set_ylabel('OD 555 [a.u.]', fontsize=6)
ax.set_title(f'slope = {popt[0]:0.5f} intercept = {popt[1]:0.2f}', fontsize=6)
plt.savefig('./viz/2024-05-28_biuret_calibration_curve.png')