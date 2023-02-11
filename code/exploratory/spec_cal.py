# %%
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/spectrophotometer_calibration.csv')

fit = np.polyfit(data['od_600nm_measured'], data['od_600nm_true'], deg=1)
fit
_data = data[data['od_600nm_measured'] > 0.45]

# plt.loglog(_data['od_600nm_measured'], (_data['od_600nm_true']-_data['od_600nm_measured']), 'o')

popt = scipy.stats.linregress(
    np.log(_data['od_600nm_measured']), np.log(_data['od_600nm_true']))
popt
od_meas = np.linspace(0.45, 1.5, 100)
fun = np.exp(popt[0] * np.log(od_meas) + popt[1])
plt.loglog(data['od_600nm_measured'], data['od_600nm_true'], 'o')
poptg
