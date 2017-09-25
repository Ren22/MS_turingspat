from sympy import *
import sympy
import numpy as np
import pandas as pd

a, b, c, d, au, u0, v0, n, M, u, v = symbols('a, b, c, d, au, u0, v0, n, M, u, v ')

f = a * (u - u0) * M ** n + b * (v - v0) - (au * (u - u0) ** 3)
g = c * (u - u0) * M ** n + d * (v - v0)

fu = diff(f, u)
fv = diff(f, v)
gu = diff(g, u)
gv = diff(g, v)

Du = 1.
Dv = 20.

l = 2 * np.pi * (fu / (2 * Du) + gv / (2 * Dv)) ** (-1 / 2.)
# print(l)

for size in np.arange(10., 120., 1.):
    l_subs = l.subs([(a, 0.2), (b, -0.4), (c, 0.6), (d, -0.8),
                     (au, 5.), (u0, 0.5), (v0, 0.5), (n, -2.), (M, 100. * (1. / size)),
                     (u, 0.5), (v, 0.5)])

    df_new = pd.DataFrame([[size, Dv / Du, l_subs]],
                          columns=['size', 'Dv/Du', 'Model #'])
    try:
        df = pd.read_csv('results/analytical_curve.csv')
        df = df.append(df_new, ignore_index=True)
    except:
        df = df_new

    df.to_csv('results/analytical_curve.csv', index=False)
    # print(result_subs)
