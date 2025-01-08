#%%
import tqdm 
import glob
import pandas as pd 
import fcsparser

# Load the file of 'valid' flow cytometry measurements and restrict to wildtype.
metadata = pd.read_csv('./flow_cytometry_sample_od600nm.csv')
metadata = metadata[metadata['strain']=='wildtype']

# Iterate through the experiment folder and load the fcs files. 
# Define constants
flow_rate = 24/60  # in µl per sec
dil_factor = 91
culture_volume = 1  # in mL

# Load the files
files = glob.glob('./202*/*.fcs')

counts = pd.DataFrame([])
for i, f in enumerate(tqdm.tqdm(files)):
        

    # Parse the metadata
    date, run_no = f.split('/')[-2].split('_')
    run_no = int(run_no[1:])
    strain, carbon, _, rep, time = f.split('/')[-1].split('_')
    rep = int(rep)
    time = float(time.split('s')[0])

    # Skip over all sample that are not wildtype
    if 'lpp' in strain:
        continue

    # Drop a spurious, large-outlying replicate for LB
    if carbon == 'LB' and  date == '2022-11-17':
        continue

    # parse the flow cytometry file and apply gating if necessary
    _, df = fcsparser.parse(f)
    n_events = len(df)

    # Determine the biomass of the particular experiment
    od = metadata[(metadata['strain'] == strain) &
                    (metadata['carbon_source'] == carbon) &
                    (metadata['run_no'] == run_no) &
                    (metadata['date'] == date)]['od_600nm'].values[0]

    # Compute the number of cells per unit biomass
    n_cells = n_events * dil_factor * \
        culture_volume * 1E3 / (flow_rate * time * od)

    # Pack the output
    _data = {'date': date,
             'run_no': run_no,
             'strain': strain,
             'carbon_source': carbon,
             'technical_rep': rep,
             'n_events': n_events,
             'cells_per_biomass': n_cells}
    _df = pd.DataFrame(_data, index=[0])
    counts = pd.concat([counts, _df], sort=False)

# %%
# Average across technical replicates
counts = counts.groupby(
    ['date', 'strain', 'carbon_source', 'run_no']).mean().reset_index()

#%%
# Load in the growth rate measurements from the mass spectrometry samples.
growth_rates = pd.read_csv('../mass_spectrometry/compiled_data/aggregated_growth_measurements.csv')
growth_rates = growth_rates[growth_rates['strain']=='wildtype']
growth_rates = growth_rates.groupby('carbon_source')['growth_rate_hr'].mean().reset_index()
 
# Create a hashmap
mapper = {g[0]:g[1] for g, _ in growth_rates.groupby(['carbon_source', 'growth_rate_hr'])}

# Add in the growth rates. 
for k, v in mapper.items():
    counts.loc[counts['carbon_source']==k, 'growth_rate_hr'] = v
counts.to_csv('../../../data/collated/experimental_cells_per_biomass.csv', index=False)





