# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()

model = cmdstanpy.CmdStanModel(
    stan_file='./cytoplasmic_density_approximation_inference.stan')

size_data = pd.read_csv(
    '../../../data/literature/collated_literature_size_data.csv')
prot_data = pd.read_csv(
    '../../../data/literature/collated_protein_per_cell.csv')
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')
dna_data = pd.read_csv(
    '../../../data/literature/collated_dna_protein_ratio.csv')

# %%
phi_cyt = ms_data[ms_data['localization'] == 'cytoplasm']
phi_rib = ms_data[ms_data['localization'] == 'ribosomal sector']

data_dict = {
    'N_size': len(size_data),
    'surface_areas': size_data['surface_area_um2'].values,
    'volume': size_data['volume_um3'].values,
    'size_lam': size_data['growth_rate_hr'].values,
    'N_dna': len(dna_data),
    'dna_protein_ratio': dna_data['DNA_protein_ratio'].values,
    'N_prot': len(prot_data),
    'prot_per_cell': prot_data['fg_protein_per_cell'].values,
    'prot_per_cell_lam': prot_data['growth_rate_hr'].values,

    'N_ms': len(phi_cyt),
    'ms_lam': phi_cyt['growth_rate_hr'].values,
    'phi_rib': phi_rib['mass_frac'].values,
    'phi_cyt': phi_cyt['mass_frac'].values
}


_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)
# prot_data = pd.read_csv('../../data/')
