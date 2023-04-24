# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import size.viz
mapper = size.viz.lit_mapper()
cor, pal = size.viz.matplotlib_style()

# Load datasets
growth_params = pd.read_csv('../../data/mcmc/growth_parameter_percentiles.csv')
model_params = pd.read_csv(
    '../../data/mcmc/perturbation_parameter_percentiles.csv')
size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
