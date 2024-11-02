"""
This script processes the mass spectrometry data collected for ppGpp perturbations
in faster growth conditions.
"""
#%%
import pandas as pd 
import tqdm

# Define the mapper given sample labels
MAPPER = {1: {"strain": "meshI", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 100,
              "replicate": 1},
          2: {"strain": "meshI", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 0,
              "replicate": 1},
          3: {"strain": "relA", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 0,
              "replicate": 1},
          4: {"strain": "relA", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 2,
              "replicate": 1},
          5: {"strain": "relA", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 4,
              "replicate": 1},
          6: {"strain": "meshI", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 100,
              "replicate": 2},
          7: {"strain": "meshI", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 0,
              "replicate": 2},
          8: {"strain": "relA", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 0,
              "replicate": 2},
          9: {"strain": "relA", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 2,
              "replicate": 2},
          10: {"strain": "relA", 
              "carbon_source": "glucoseCAA",
              "date": "2024-08-27",
              "inducer_conc": 4,
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
standard = standard[(standard['growth_rate_hr'] >= 0.75) & 
                    (standard['growth_rate_hr'] <= 0.95)]
standard = standard.groupby('gene_name')['mass_frac'].mean().reset_index()
std_map = standard.set_index(['gene_name'])['mass_frac'].to_dict()


df = pd.DataFrame([])
unmapped = []
raw_df = pd.read_csv(f'./raw/all_replicate_intensities.csv')
std_norm = raw_df['Normalized_proportion_1']
for i in range(len(MAPPER)):
    raw_df[f'mass_frac_{i+1}'] = raw_df['Intensity'].values * raw_df[f'Normalized_proportion_{i+1}'].values / raw_df[f'Normalized_proportion_{i+1}'].sum()
    raw_df[f'mass_frac_{i+1}'] /= raw_df[f'mass_frac_{i+1}'].sum()

for g, d in tqdm.tqdm(raw_df.groupby('Protein')):  
    for i in range(len(MAPPER)):
        _mapper = MAPPER[i+1]
        _df_dict = {k:v for k, v in _mapper.items()}
        gene_names = uniprot_mapper[g]['gene_names'].split()
        _df_dict['entry'] = g
        _df_dict['name'] = gene_names[0] 
        _df_dict['synonyms'] = ' '.join(gene_names[1:])
        _df_dict['intensity'] = d[f'Intensity_{i+1}'].values[0]
        _df_dict['norm_intensity'] = d[f'Normalized_proportion_{i+1}'].values[0]
        _df_dict['common_intensity'] = d['Intensity'].values[0]
        _df_dict['mass_frac'] = d[f'mass_frac_{i+1}'].values[0]
        df = pd.concat([df, pd.DataFrame(_df_dict, index=[0])])
df.to_csv('../processed_mass_fractions_ppGpp_glucoseCAA.csv', index=False)

