# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import size.viz
cor, pal = size.viz.matplotlib_style()

voldata = pd.read_csv(
    '../../data/mcmc/wildtype_hyperparameter_size_summary.csv')
grdata = pd.read_csv('../../data/mcmc/wildtype_growth_rate_summary.csv')


def width_pred(lam,
               delta=0.025,
               k=1,
               alpha=4,
               kappa=0.04,
               phi_0=0.15):
    return (12 * alpha * k * delta * (1 - phi_0 + kappa * lam))/((phi_0 - kappa * lam) * (3 * alpha - 1))


# %%
lam_range = np.linspace(0, 1.99, 200)
w = width_pred(lam_range)
plt.plot(lam_range, w, '--')

for g, d in grdata.groupby(['carbon_source']):
    vol = voldata[(voldata['carbon_source'] == g) &
                  (voldata['parameter'] == 'width_um')]
    plt.plot(d['median'], vol['median'], 'o')

# %%
voldata
