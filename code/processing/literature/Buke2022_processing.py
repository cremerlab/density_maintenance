# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import glob

files = glob.glob('../../../data/literature/Buke2022/*.mat')
mat = scipy.io.loadmat(files[0], squeeze_me=True)
mat
# %%
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')


outer = ms_data[ms_data['localization'] == 'outer membrane']
inner = ms_data[ms_data['localization'] == 'inner membrane']

outer['rho'] = outer['mass_fg'] / outer['surface_area']
inner['rho'] = inner['mass_fg'] / inner['surface_area']


plt.plot(outer['growth_rate_hr'], outer['rho'], 'o', color='dodgerblue')
plt.plot(inner['growth_rate_hr'], inner['rho'], 'o', color='rebeccapurple')
plt.ylim([0, 6])
