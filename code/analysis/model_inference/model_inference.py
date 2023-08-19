# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()

size_data = pd.read_csv(
    '../../data/literature/collated_literature_size_data.csv')
size_data = size_data[size_data['source'] != 'Si et al. 2017']

prot_data = pd.read_csv('../../data/literature/collated_protein_per_cell.csv')

ms_data = pd.read_csv(
    '../../data/literature/collated_mass_fractions_empirics.csv')
mem_data = ms_data[ms_data['localization'] == 'membrane']
peri_data = ms_data[ms_data['localization'] == 'periplasm']
phiRb_data = pd.read_csv(
    '../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')
phiRb_data = phiRb_data[phiRb_data['source'] != 'Si et al. 2017']

# %%
model = cmdstanpy.CmdStanModel(stan_file='mass_spec_empirics_uncertainty.stan')
