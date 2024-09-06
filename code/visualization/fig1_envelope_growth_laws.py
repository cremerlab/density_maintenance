#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

# Load our data 
data = pd.read_csv('../processing/mass_spectrometry/mass_spec_sector_assignments.csv')
data = data[data['strain']=='wildtype']

# Load the literature data 
lit_data = pd.read_csv('../../')
data.head()
