#%%

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
ms_conf = pd.read_csv(
    '../../data/mass_spectrometry/periplasm_fractionation_confirmation.csv')

grouped = ms_conf[ms_conf['localization']=='periplasm'][['fractional_peri', 'fractional_cyto']].sum().reset_index()
grouped

fig, ax = plt.subplots(1,1, figsize=(1, 1))
ax.bar(['supernatant', 'pellet'], grouped[0].values)
ax.set_ylabel('fraction of total intensity', fontsize=6)