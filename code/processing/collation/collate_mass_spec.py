#%%
import pandas as pd 
################################################################################
## Merge individual gene copy numbers, only for wildtype
################################################################################

# Load the full literature data 
lit_data = pd.read_csv('../literature/compiled_data/compiled_literature_mass_fractions.csv')

# Load experimental data
exp_data = pd.read_csv('../../../data/collated/experimental_mass_spectrometry.csv')

# Unify column names
exp_data = exp_data[exp_data['strain']=='wildtype']
exp_data['source'] = 'This Study'
exp_data.drop(columns=['date', 'inducer_conc', 'strain'], inplace=True)
lit_data['replicate'] = 0
lit_data.rename(columns={'condition':'carbon_source',
                         'cog_class':'cog_category'}, inplace=True)
merged = pd.concat([exp_data, lit_data])
merged.to_csv('../../../data/collated/merged_mass_spectrometry.csv', index=False)


