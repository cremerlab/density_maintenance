# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Si et al. 2017']

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
size_emp = pd.read_csv(
    '../../data/mcmc/size_data_empirical_summaries_wide.csv')
ms_emp = pd.read_csv(
    '../../data/mcmc/mass_spec_empirical_summaries_wide.csv')
theo_lam = pd.read_csv(
    '../../data/mcmc/theory_growth_rate_prediction_summaries.csv')
theo_phiRb = pd.read_csv(
    '../../data/mcmc/theory_phiRb_prediction_summaries.csv')

# %%
