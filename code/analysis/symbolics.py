# %%
import sympy as sp
k = sp.Symbol('k')
Lam = sp.Symbol('{{\Lambda}}')
w_min = sp.Symbol('{{w_{min}}}')
phi = sp.Symbol('\phi')
phi_max = 1/(1/k * (w_min/Lam - 1) + 1)
dphi = phi - phi_max
taylor = w_min - (Lam * k * dphi / phi_max**2) * (dphi/phi_max - 1)
taylor.simplify()
