#%%
import numpy as np 
import pandas as pd 

# Load the cell count data
flow_data = pd.read_csv('../../../data/summaries/summarized_cell_counts.csv')
flow_data = flow_data[flow_data['strain']=='wildtype']
flow_data_agg = flow_data.groupby(['strain', 'carbon_source'])['cells_per_biomass'].agg((np.mean, np.std)).reset_index()
flow_data_agg.rename(columns={'mean': 'mean_cells_per_biomass', 'std': 'std_cells_per_biomass'}, inplace=True)

# Load the total protein data
prot_data = pd.read_csv('../../../data/summaries/summarized_total_protein.csv')
prot_data = prot_data[prot_data['overexpression']=='none']
prot_data_agg = prot_data.groupby(['carbon_source'])['ug_prot_per_biomass'].agg((np.mean, np.std)).reset_index()
prot_data_agg.rename(columns={'mean': 'mean_ug_prot_per_biomass', 'std': 'std_ug_prot_per_biomass'}, inplace=True)

# Load the growth rate data collected by roshali
gr_data = pd.read_csv('../../../data/summaries/summarized_growth_measurements.csv')
gr_data = gr_data[(gr_data['strain']=='wildtype') & (gr_data['overexpression']=='none') &
                  (gr_data['temperature']==37)]
gr_data_agg = gr_data.groupby(['carbon_source'])['growth_rate_hr'].agg((np.mean, np.std)).reset_index()
gr_data_agg.rename(columns={'mean': 'mean_growth_rate_hr', 'std': 'std_growth_rate_hr'}, inplace=True)

# Merge the data
merged_data = pd.merge(prot_data_agg, flow_data_agg, on='carbon_source', how='inner')
merged_data = pd.merge(merged_data, gr_data_agg, on='carbon_source', how='inner')
merged_data['fg_prot_per_cell'] = 1E9 * merged_data['mean_ug_prot_per_biomass'] / merged_data['mean_cells_per_biomass']

# Assume uncorrelated errors and propagate
merged_data['err_fg_prot_per_cell'] =  merged_data['fg_prot_per_cell'] *\
    np.sqrt( (merged_data['std_ug_prot_per_biomass'] / merged_data['mean_ug_prot_per_biomass'])**2 +\
             (merged_data['std_cells_per_biomass'] / merged_data['mean_cells_per_biomass'])**2 )
merged_data = merged_data[['carbon_source', 'fg_prot_per_cell', 'err_fg_prot_per_cell', 'mean_growth_rate_hr', 'std_growth_rate_hr']]
merged_data.to_csv('../../../data/bulk_protein_per_cell.csv', index=False)
