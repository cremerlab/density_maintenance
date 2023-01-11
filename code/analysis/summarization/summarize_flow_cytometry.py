# %%
import numpy as np
import pandas as pd
import tqdm
import fcsparser
import glob

# Define constants
flow_rate = 24/60  # in µl per sec
dil_factor = 91
culture_volume = 1  # in mL

# Load the file containing the info for the OD harvest
sample_ods = pd.read_csv(
    '../../../data/flow_cytometry/flow_cytometry_sample_od600nm.csv')

# Load the files
files = glob.glob('../../../data/flow_cytometry/2022*/*.fcs')

counts = pd.DataFrame([])
for i, f in enumerate(tqdm.tqdm(files)):

    # Parse the metadata
    date = f.split('/')[-2]
    strain, carbon, _, rep, time = f.split('/')[-1].split('_')
    rep = int(rep)
    time = float(time.split('s')[0])

    # parse the flow cytometry file and apply gating if necessary
    _, df = fcsparser.parse(f)
    n_events = len(df)

    # Determine the biomass of the particular experiment
    od = sample_ods[(sample_ods['strain'] == strain) &
                    (sample_ods['carbon_source'] == carbon) &
                    (sample_ods['date'] == date)]['od_600nm'].values[0]

    # Compute the number of cells per unit biomass
    n_cells = n_events * dil_factor * \
        culture_volume * 1E3 / (flow_rate * time * od)

    # Pack the output
    _data = {'strain': strain,
             'carbon_source': carbon,
             'overexpression': 'none',
             'inducer': 'none',
             'inducer_conc': 0,
             'technical_rep': rep,
             'n_events': n_events,
             'gated': False,
             'cells_per_biomass': n_cells}
    _df = pd.DataFrame(_data, index=[0])
    counts = pd.concat([counts, _df], sort=False)

# %%
counts.to_csv('../../../data/summaries/flow_cytometry_counts.csv', index=False)
