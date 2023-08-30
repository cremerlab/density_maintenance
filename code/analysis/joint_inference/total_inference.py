# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
import arviz as az
import cmdstanpy
cor, pal = size.viz.matplotlib_style()

model = cmdstanpy.CmdStanModel(stan_file='./total_inference.stan')

# Load the literature datasets
lit_size = pd.read_csv(
    '../../../data/literature/collated_literature_size_data.csv')
lit_prot = pd.read_csv(
    '../../../data/literature/collated_protein_per_cell.csv')
lit_drymass = pd.read_csv(
    '../../../data/literature/collated_drymass_densities.csv')
ms_data = pd.read_csv(
    '../../../data/literature/collated_mass_fractions_empirics.csv')
mem = ms_data[ms_data['localization'] == 'membrane']
peri = ms_data[ms_data['localization'] == 'periplasm']

lit_phiRb = pd.read_csv(
    '../../../data/literature/Chure2023/chure2023_collated_mass_fractions.csv')

# Load our wildtype data
total_protein = pd.read_csv(
    '../../../data/summaries/summarized_total_protein.csv')
total_protein = total_protein[(total_protein['strain'] == 'wildtype') & (
    total_protein['overexpression'] == 'none')]

total_peri = pd.read_csv(
    '../../../data/summaries/summarized_periplasmic_protein.csv')
total_peri = total_peri[(total_peri['strain'] == 'wildtype') & (
    total_peri['overexpression'] == 'none')]

total_mem = pd.read_csv(
    '../../../data/summaries/summarized_membrane_protein.csv')
total_mem = total_mem[(total_mem['strain'] == 'wildtype')
                      & (total_mem['overexpression'] == 'none')]

total_rna = pd.read_csv('../../../data/summaries/summarized_total_rna.csv')
total_rna = total_rna[(total_rna['strain'] == 'wildtype')
                      & (total_rna['overexpression'] == 'none')]

growth_data = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')
growth_data = growth_data[(growth_data['strain'] == 'wildtype') &
                          (growth_data['overexpression'] == 'none')]

size_data = pd.read_csv(
    '../../../data/summaries/summarized_size_measurements.csv')
size_data = size_data[(size_data['strain'] == 'wildtype') &
                      (size_data['overexpression'] == 'none')]

flow_data = pd.read_csv('../../../data/summaries/summarized_cell_counts.csv')
flow_data = flow_data[flow_data['cells_per_biomass'] > 1E6]

# %%
# Apply a condition mapper
mapper = {'acetate': 1,
          'sorbitol': 2,
          'glycerol': 3,
          'glucose':  4,
          'glucoseCAA': 5,
          'LB': 6}

for d in [total_protein, total_mem, total_rna, total_peri, growth_data, size_data, flow_data]:
    for k, v in mapper.items():
        d.loc[d['carbon_source'] == k, 'idx'] = v
    d.dropna(inplace=True)
    d['idx'] = d['idx'].values.astype(int)

# %%
N_pred = 500
pred_phiRb = np.linspace(0.06, 0.35, N_pred)


# Set up the data dictionary
data_dict = {

    'N_pred': N_pred,
    'pred_phiRb':  pred_phiRb,

    'N_size_lit': len(lit_size),
    'lit_size_lam': lit_size['growth_rate_hr'].values,
    'lit_SAV': lit_size['surface_to_volume'].values,

    'N_prot_lit': len(lit_prot),
    'lit_prot_lam': lit_prot['growth_rate_hr'].values,
    'lit_prot_per_cell': lit_prot['fg_protein_per_cell'].values,

    'N_drymass_lit': len(lit_drymass),
    'lit_drymass_lam': lit_drymass['growth_rate_hr'].values,
    'lit_drymass_density': lit_drymass['drymass_density_fg_fL'].values,

    'N_ms': len(mem),
    'lit_ms_lam': mem['growth_rate_hr'].values,
    'lit_phi_mem': mem['mass_frac'].values,
    'lit_phi_peri': peri['mass_frac'].values,

    'N_phiRb_lit': len(lit_phiRb),
    'lit_phiRb_lam': lit_phiRb['growth_rate_hr'].values,
    'lit_phiRb': lit_phiRb['mass_fraction'].values,

    'N_prot': len(total_protein),
    'J_prot': total_protein['idx'].max(),
    'prot_idx': total_protein['idx'].values,
    'prot_per_biomass': total_protein['ug_prot_per_biomass'].values,

    'N_peri': len(total_peri),
    'J_peri': total_peri['idx'].max(),
    'peri_idx': total_peri['idx'].values,
    'peri_per_biomass': total_peri['ug_prot_per_biomass'].values,

    'N_mem': len(total_mem),
    'J_mem': total_mem['idx'].max(),
    'mem_idx': total_mem['idx'].values,
    'mem_per_biomass': total_mem['ug_prot_per_biomass'].values,

    'N_rna': len(total_rna),
    'J_rna': total_rna['idx'].max(),
    'rna_idx': total_rna['idx'].values,
    'rna_per_biomass': total_rna['ug_rna_per_biomass'].values,

    'N_growth': len(growth_data),
    'J_growth': growth_data['idx'].max(),
    'growth_idx': growth_data['idx'].values,
    'growth_rates': growth_data['growth_rate_hr'].values,

    'N_size': len(size_data),
    'J_size': size_data['idx'].max(),
    'size_idx': size_data['idx'].values,
    'surface_areas': size_data['surface_area'],
    'volumes': size_data['volume'],

    'N_flow': len(flow_data),
    'J_flow': flow_data['idx'].max(),
    'flow_idx': flow_data['idx'].values,
    'cells_per_biomass':  flow_data['cells_per_biomass'].values
}

# %%
_samples = model.sample(data=data_dict)
samples = az.from_cmdstanpy(_samples)

# %%
pred = samples.posterior.pred_SAV.to_dataframe().reset_index()
lower, upper = [], []
for g, d in pred.groupby('pred_SAV_dim_0'):
    perc = np.percentile(d['pred_SAV'], (2.5, 97.5))
    lower.append(perc[0])
    upper.append(perc[1])

fig, ax = plt.subplots(1, 1)
ax.fill_between(pred_phiRb, lower, upper, color=cor['light_blue'], alpha=0.5)

for g, d in samples.posterior.phiRb_SAV.to_dataframe().reset_index().groupby('phiRb_SAV_dim_0'):
    med_val = np.median(d['phiRb_SAV'])
    perc = np.percentile(d['phiRb_SAV'].values, (2.5, 97.5))
    phiRb_val = lit_phiRb['mass_fraction'].values[g]
    fmt = size.viz.style_point(lit_phiRb['source'].values[g])
    ax.vlines(phiRb_val, *perc, linewidth=1,
              color=cor['primary_black'], alpha=0.25)
    ax.plot(phiRb_val, med_val, **fmt)

# %%
plt.plot(lit_size['growth_rate_hr'], lit_size['surface_area_um2'], 'o')
plt.ylim([0, 10])
