# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('../../data/literature/collated_mass_fractions.csv')
abcs = data[data['go_terms'].str.contains('GO:0055052')]
abcs
