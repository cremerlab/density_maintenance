#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load the many literature data sets
lit_ms = pd.read_csv('../../data/collated/compiled_literature_allocation_assignments_wide.csv')
lit_rp = pd.read_csv('../../data/collated/collated_literature_rna_to_protein.csv')
lit_size = pd.read_csv('../../data/collated/collated_literature_size_data.csv')
lit_prot = pd.read_csv('../../data/collated/collated_literature_total_protein.csv')
lit_rna = pd.read_csv('../../data/collated/collated_literature_total_RNA.csv')

# Load our data sets
our_ms = pd.read_csv('../../data/collated/experimental_mass_spectrometry.csv')
our_size = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
our_rp = pd.read_csv('../../data/collated/experimental_rna_protein_per_cell.csv')

