from sympy import *
import sympy
import numpy as np
import pandas as pd

### Used instability criteria ###
#       0) (sqrt(fu * gv - fv * gu) + sqrt(-fv * gu)) ** 2)/fu ** 2
#       1) fu+gv < 0 = {0},
#       2) fu*gv-fv*gu > 0 = {1},
#       3) dfu+gv > 0 = {2},
#       4)(dfu+gv)^2-4d(fugv-fvgu) > 0 = {3}

a, b, c, d, au, u0, v0, n, M, u, v = symbols('a, b, c, d, au, u0, v0, n, M, u, v ')

f = a * (u - u0) * M ** n + b * (v - v0) - (au * (u - u0) ** 3)
g = c * (u - u0) + d * (v - v0)
fu = diff(f, u)
fv = diff(f, v)
gu = diff(g, u)
gv = diff(g, v)

Du = 1.
Dv = 20.

numerator = (sympy.sqrt(fu * gv - fv * gu) + sympy.sqrt(-fv * gu)) ** 2
denominator = (fu ** 2)

criter0 = numerator / denominator
criter1 = fu + gv
criter2 = fu * gv - fv * gu
criter3 = (Dv / Du) * fu + gv
criter4 = ((Dv / Du) * fu + gv) ** 2 - 4 * (Dv / Du) * (fu * gv - fv * gu)

# print(criter0)
for size in np.arange(10., 200., 1.):
    params_list = [(a, 0.65), (b, -0.95), (c, 0.6), (d, -0.8),
                   (au, 5.), (u0, 0.5), (v0, 0.5), (n, -2.), (M, 100. * (1. / size)),
                   (u, 0.5), (v, 0.5)]
    criter0_subs, criter1_subs, criter2_subs, criter3_subs, criter4_subs = \
        criter0.subs(params_list), \
        criter1.subs(params_list), \
        criter2.subs(params_list), \
        criter3.subs(params_list), \
        criter4.subs(params_list)

    df_new = pd.DataFrame([[size, Dv / Du, criter0_subs, criter1_subs, criter2_subs, criter3_subs, criter4_subs]],
                          columns=['size', 'Dv/Du', 'Cr.0: [val]<Dv/Du', 'Cr.1:[val]<0',
                                   'Cr.2:[val]>0', 'Cr.3:[val]>0', 'Cr.4:[val]>0'])
    try:
        df = pd.read_csv('results/all_instability_crits.csv')
        df = df.append(df_new, ignore_index=True)
    except:
        df = df_new

    df.to_csv('results/all_instability_crits.csv', index=False)
    # print(result_subs)
