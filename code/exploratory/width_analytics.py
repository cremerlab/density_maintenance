# %%
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import size.viz
import palettable.cmocean as cmocean
cor, pal = size.viz.matplotlib_style()
cmap = cmocean.sequential.Ice_5.mpl_colors
phi_range = np.linspace(0.001, 1, 200)
k_vals = [0.01, 0.1, 0.5, 1]
delta = 0.020
alpha = 3


fig, ax = plt.subplots(1, 2, figsize=(5, 2))
for i, k in enumerate(k_vals):
    w = (12 * alpha * delta / (3 * alpha - 1)) * \
        (1 + k * (1 - phi_range)/phi_range)
    ax[0].plot(phi_range, w, label=k, lw=1, color=cmap[i])
ax[0].legend(title='k')
ax[0].set_ylim([0, 1.5])
ax[0].set_ylabel('$w$\ncell width [µm]')
ax[0].set_xlabel('periplasmic mass fraction\n'+r'$\varphi_M$')

alpha = 3
delta = 0.024
Lam = 12 * alpha * delta / (3 * alpha - 1)
w_min = 0.5
w_max = 1
delta_w = w_max - w_min

k_vals = [0.01, 0.1, 0.5, 1, 1.5]
for k in k_vals:
    phi_min = Lam * k / (w_max + Lam * (k - 1))
    phi_max = Lam * k / (w_min + Lam * (k - 1))
    phi_range = np.linspace(phi_min, phi_max, 100)
    slope = Lam * k / (phi_max * phi_min)
    pred = w_max - slope * (phi_range - phi_min)

    # pred = w_max - (Lam * k / 0.06) * (-1 + phi_range / 0.005)
    pred = w_max + (w_min + Lam * (k - 1)) * \
        (1 - phi_range * (w_max + Lam * (k - 1))/(Lam * k))

    ax[1].plot(phi_range, pred)

# %%
k = sp.Symbol('k')
delta = sp.Symbol('\delta')
alpha = sp.Symbol(r'\alpha')
wmax = sp.Symbol('{{w_{max}}}')
wmin = sp.Symbol('{{w_{min}}}')
phi = sp.Symbol('\phi')
theta = sp.Symbol(r'\theta')
# theta = 12 * k * delta * alpha / (3 * alpha - 1)
# wmin = 1/2
# wmax = 1
phi_min = theta / (wmax + theta * (1 - 1/k))
phi_max = theta / (wmin + theta * (1 - 1/k))
eq = wmax + (theta / phi_max) * (1 - phi / phi_min)
eq.to_latex()
