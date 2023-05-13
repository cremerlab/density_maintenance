# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/literature/collated_mass_fractions.csv')
# %%
# data = data[(data['dataset_name'] == 'Schmidt et al. 2016')
# & (data['condition'] == 'lb_miller')]

cytoplasm = data[(data['envelope'] == False) & (data['cog_class'].isin(
    ['information storage and processing', 'metabolism']))]
envelope = data[(data['envelope'] == True)]

fig, ax = plt.subplots(1, 4)
