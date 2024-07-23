"""
This script processes the mass spectrometry data collected for all 32 samples.
"""
#%%
import numpy as np 
import pandas as pd 
import tqdm

# Define the mapper given sample labels
MAPPER_REP1 = {1: {"strain": "wildtype", 
              "carbon_source": "glucose",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 1},
          2: {"strain": "wildtype", 
              "carbon_source": "glucose+acetate",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 1},
          3: {"strain": "wildtype", 
              "carbon_source": "glycerol",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 1},
          4: {"strain": "wildtype", 
              "carbon_source": "sorbitol",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 1},
          5: {"strain": "wildtype", 
              "carbon_source": "glucoseCAA",
              "date": "2024-05-24",
              "inducer_conc": 0,
              "replicate": 1},
          6: {"strain": "wildtype", 
              "carbon_source": "LB",
              "date": "2024-05-24",
              "inducer_conc": 0,
              "replicate": 1},
          7: {"strain": "wildtype", 
              "carbon_source": "acetate",
              "date": "2024-05-25",
              "inducer_conc": 0,
              "replicate": 1},
          8: {"strain": "meshI", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 0,
              "replicate": 1},
          9: {"strain": "meshI", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 100,
              "replicate": 1},
          10: {"strain": "relA", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 0,
              "replicate": 1},
          11: {"strain": "relA", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 2,
              "replicate": 1},
          12: {"strain": "relA", 
              "carbon_source": "glucose",
              "date": "2024-05-25",
              "inducer_conc": 4,
              "replicate": 1},
          13: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-25",
              "inducer_conc": 0,
              "replicate": 1},
          14: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-28",
              "inducer_conc": 1,
              "replicate": 1},
          15: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-29",
              "inducer_conc": 3,
              "replicate": 1},
          16: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-25",
              "inducer_conc": 5,
              "replicate": 1},
        }
MAPPER_REP2 = {1: {"strain": "wildtype", 
              "carbon_source": "glucose",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 2},
          2: {"strain": "wildtype", 
              "carbon_source": "glucose+acetate",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 2},
          3: {"strain": "wildtype", 
              "carbon_source": "glycerol",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 2},
          4: {"strain": "wildtype", 
              "carbon_source": "sorbitol",
              "date": "2024-05-22",
              "inducer_conc": 0,
              "replicate": 2},
          5: {"strain": "wildtype", 
              "carbon_source": "glucoseCAA",
              "date": "2024-05-24",
              "inducer_conc": 0,
              "replicate": 2},
          6: {"strain": "wildtype", 
              "carbon_source": "LB",
              "date": "2024-05-24",
              "inducer_conc": 0,
              "replicate": 2},
          7: {"strain": "wildtype", 
              "carbon_source": "acetate",
              "date": "2024-05-25",
              "inducer_conc": 0,
              "replicate": 2},
          8: {"strain": "meshI", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 0,
              "replicate": 2},
          9: {"strain": "meshI", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 100,
              "replicate": 2},
          10: {"strain": "relA", 
              "carbon_source": "glucose",
              "date": "2024-05-23",
              "inducer_conc": 0,
              "replicate": 2},
          11: {"strain": "relA", 
              "carbon_source": "glucose",
              "date": "2024-05-26",
              "inducer_conc": 2,
              "replicate": 2},
          12: {"strain": "relA", 
              "carbon_source": "glucose",
              "date": "2024-05-25",
              "inducer_conc": 4,
              "replicate": 2},
          13: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-25",
              "inducer_conc": 0,
              "replicate": 2},
          14: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-28",
              "inducer_conc": 1,
              "replicate": 2},
          15: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-26",
              "inducer_conc": 3,
              "replicate": 2},
          16: {"strain": "lacZ", 
              "carbon_source": "glucose",
              "date": "2024-05-25",
              "inducer_conc": 5,
              "replicate": 2},
        }
uniprot = pd.read_csv('./raw/uniprot_identifiers.csv')
uniprot_mapper = {}
for i in range(len(uniprot)):
    g = uniprot.iloc[i]
    uniprot_mapper[g['Entry name']] = {'gene_names': g['Gene names']} 

#%%
# Load the two files
files = ['replicate1_intensities.csv',
         'replicate2_intensities.csv']

standard = pd.read_csv('../../../../data/literature/compiled_mass_fractions.csv')
standard = standard[(standard['strain']=='NCM3722') & 
                    (standard['growth_rate_hr'] >= 0.90) & 
                    (standard['growth_rate_hr'] <= 1.2)]
standard = standard.groupby('gene_name')['mass_frac'].mean().reset_index()
std_map = standard.set_index(['gene_name'])['mass_frac'].to_dict()

#%%
df = pd.DataFrame([])
unmapped = []
for mapper, file in zip([MAPPER_REP1, MAPPER_REP2], files):
    raw_df = pd.read_csv(f'./raw/{file}')
    std_norm = raw_df['Normalized_proportion_1']
    for j in range(16):
        raw_df[f'mass_frac_{j+1}'] = raw_df['Intensity'].values * raw_df[f'Normalized_proportion_{j+1}'].values / raw_df[f'Normalized_proportion_{j+1}'].sum()
        raw_df[f'mass_frac_{j+1}'] /= raw_df[f'mass_frac_{j+1}'].sum()
        raw_df[f'relative_norm_prop_{j+1}'] = raw_df[f'Normalized_proportion_{j+1}'].values / std_norm.values
    for g, d in tqdm.tqdm(raw_df.groupby('Protein')):  
        for j in range(16):
            _mapper = mapper[j+1]
            _df_dict = {k:v for k, v in _mapper.items()}
            gene_names = uniprot_mapper[g]['gene_names'].split()
            try:
                _std = std_map[gene_names[0]]
            except KeyError:
                unmapped.append(gene_names[0])
                _std = 0
            _df_dict['entry'] = g
            _df_dict['name'] = gene_names[0] 
            _df_dict['synonyms'] = ' '.join(gene_names[1:])
            _df_dict['intensity'] = d[f'Intensity_{j+1}'].values[0]
            _df_dict['norm_intensity'] = d[f'Normalized_proportion_{j+1}'].values[0]
            _df_dict['common_intensity'] = d['Intensity'].values[0]
            _df_dict['mass_frac'] = d[f'mass_frac_{j+1}'].values[0]
            _df_dict['rel_mass_frac'] = _std * d[f'relative_norm_prop_{j+1}'].values[0]
            df = pd.concat([df, pd.DataFrame(_df_dict, index=[0])])
df.to_csv('../processed_mass_fractions.csv', index=False)

