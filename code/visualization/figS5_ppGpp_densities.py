#%%
import pandas as pd
import matplotlib.pyplot as plt
import size.viz 
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('../../data/collated/aggregated_experimental_data.csv')
data.head()