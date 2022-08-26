# %%
from sqlite3 import adapt
import numpy as np
import pandas as pd
import cmdstanpy
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/wt_cell_sizes.csv')
model = cmdstanpy.CmdStanModel(stan_file='stan_models/size_inference.stan')
glu = data[data['carbon_source'] == 'glucose']
glu['cell_id'] += 1
# %%
J = 1
N = int(len(glu))
widths = glu['width_mean'].values.astype(float)
lengths = glu['length'].values.astype(float)

# %%
data_dict = {
    'N': N,
    'widths': widths,
    'lengths': lengths}
samples = model.sample(data_dict)  # , adapt_delta=0.99)
print(samples.diagnose())
# %%
samples.posterior.to_dataframe()
# %%
samples.summary()

# %%
