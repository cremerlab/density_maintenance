# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmdstanpy
import arviz as az
import size.viz
cor, pal = size.viz.matplotlib_style()
model = cmdstanpy.CmdStanModel(stan_file='./data_inference.stan')

# Load the datasets
prot = pd.read_csv('../../../data/summaries/summarized_total_protein.csv')
mem = pd.read_csv('../../../data/summaries/summarized_membrane_protein.csv')
peri = pd.read_csv(
    '../../../data/summaries/summarized_periplasmic_protein.csv')
rna = pd.read_csv('../../../data/summaries/summarized_total_rna.csv')
sizes = pd.read_csv('../../../data/summaries/summarized_size_measurements.csv')
flow = pd.read_csv('../../../data/summaries/summarized_cell_counts.csv')
flow = flow[flow['cells_per_biomass'].values > 1E6]
growth = pd.read_csv(
    '../../../data/summaries/summarized_growth_measurements.csv')

# %%
# Restrict to wildtype
prot = prot[(prot['strain'] == 'wildtype') &
            (prot['overexpression'] == 'none')]
mem = mem[(mem['strain'] == 'wildtype') & (mem['overexpression'] == 'none')]
peri = peri[(peri['strain'] == 'wildtype') &
            (peri['overexpression'] == 'none')]
rna = rna[(rna['strain'] == 'wildtype') & (rna['overexpression'] == 'none')]
sizes = sizes[(sizes['strain'] == 'wildtype') &
              (sizes['overexpression'] == 'none')]
sizes = sizes[sizes['carbon_source'] != 'ezMOPS']
growth = growth[(growth['strain'] == 'wildtype') &
                (growth['overexpression'] == 'none')]

# Add a condition mapper
mapper = {'acetate': 1, 'sorbitol': 2, 'glycerol': 3,
          'glucose': 4, 'glucoseCAA': 5, 'LB': 6}

for d in [prot, mem, peri, rna, sizes, growth, flow]:
    d['idx'] = [mapper[k] for k in d['carbon_source'].values]
    d['idx'] = d['idx'].values.astype(int)
# %%
# Sample the model
data_dict = {
    'N_prot': len(prot),
    'J_prot': prot['idx'].max(),
    'prot_idx': prot['idx'].values,
    'prot_per_biomass': prot['ug_prot_per_biomass'].values,

    'N_peri': len(peri),
    'J_peri': peri['idx'].max(),
    'peri_idx': peri['idx'].values,
    'peri_per_biomass': peri['ug_prot_per_biomass'].values,

    'N_mem': len(mem),
    'J_mem': mem['idx'].max(),
    'mem_idx': mem['idx'].values,
    'mem_per_biomass': mem['ug_prot_per_biomass'].values,

    'N_rna': len(rna),
    'J_rna': rna['idx'].max(),
    'rna_idx': rna['idx'].values,
    'rna_per_biomass': rna['ug_rna_per_biomass'].values,

    'N_growth': len(growth),
    'J_growth': growth['idx'].max(),
    'growth_idx': growth['idx'].values,
    'growth_rates': growth['growth_rate_hr'].values,

    'N_size': len(sizes),
    'J_size': sizes['idx'].max(),
    'size_idx': sizes['idx'].values,
    'widths': sizes['width_median'].values,
    'lengths': sizes['length'].values,
    'volumes': sizes['volume'].values,
    'surface_areas': sizes['surface_area'].values,
    'aspect_ratios': sizes['aspect_ratio'].values,

    'N_flow': len(flow),
    'J_flow': flow['idx'].max(),
    'flow_idx': flow['idx'].values,
    'cells_per_biomass': flow['cells_per_biomass'].values,
}

_samples = model.sample(data=data_dict, adapt_delta=0.99)
samples = az.from_cmdstanpy(_samples)
# %%
samples.posterior
