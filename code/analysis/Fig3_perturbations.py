# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.sequential import Blues_7
import size.viz
import size.analytical
cor, pal = size.viz.matplotlib_style()

# Load datasets
params = pd.read_csv('../../data/mcmc/parameter_percentiles.csv')
singulars = pd.read_csv('../../data/mcmc/singular_parameter_percentiles.csv')
shapes = pd.read_csv('../../data/mcmc/shape_posterior_kde.csv')
posts = pd.read_csv('../../data/mcmc/model_posterior_kde.csv')
pred = pd.read_csv('../../data/mcmc/predicted_scaling_lppwt.csv')

# Restrict to wildtype
params = params[~((params['overexpression']!='none') * (params['overexpression']==0))]
posts = posts[~((posts['overexpression']!='none') * (posts['overexpression']==0))]         
posts.rename(columns={'inducer_conc_ng_mL': 'inducer_conc'})
shapes = shapes[~((shapes['overexpression']!='none') * (shapes['overexpression']==0))]         
pred_err = pred[pred['interval'] != 'median']
posts = pd.concat([posts, shapes], sort=False)
singular_medians = singulars[singulars['interval'] == 'median']
singular_errs = singulars[singulars['interval'] != 'median']
medians = params[params['interval'] == 'median']
errs = params[params['interval'] != 'median']