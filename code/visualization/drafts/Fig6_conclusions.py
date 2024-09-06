# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

# size_data = pd.read_csv('../../data/literature/collated_literature_size_data.csv')
# fig, ax = plt.subplots(1, 1)
# ax.plot(size_data['width_um'], size_data['length_um'], '.')
kdes = pd.read_csv('../../data/mcmc/literature_model_params_kde.csv')
kdes = kdes[(kdes['volume_scale'] == 'linear_width') &
            (kdes['model'] == 'const_phi_mem')]
