import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from sympy import *
import sympy
import pandas as pd

a, b, c, d, au, u0, v0, n, M, u, v = symbols('a, b, c, d, au, u0, v0, n, M, u, v')

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

# Parameters range definition
a_range = np.arange(.1, .2, .1)
b_range = np.arange(-.9, -1., -.1)
d_range = np.arange(-.8, -1., -.1)
sys_size = 50.

# Matrix definition
Z = np.zeros((len(a_range), len(b_range), len(d_range)))
# Z[0, 1:] = [np.round(elem, 2) for elem in b_range]  # all b_range assigned to 0 row
# Z[1:, 0] = [np.round(elem, 2) for elem in a_range]  # all a_range assigned to 0 column

ja = 0;
jb = 0;
jd = 0;

for ia in a_range:
    for ib in b_range:
        for id in d_range:
            if (ia == 0): ia = 1e-10
            if (ib == 0): ib = 1e-10
            if (id == 0): id = 1e-10
            res = criter0.subs([(a, ia), (b, ib), (c, 0.6), (d, id),
                                (au, 5.), (u0, 0.5), (v0, 0.5), (n, -2.), (M, 100. * (1. / sys_size)),
                                (u, 0.5), (v, 0.5)])
            if res.is_Mul:
                res = 2 * Dv / Du
            Z[ja, jb, jd] = res
            print Z
            # print Z[ja, jb, jd]
            jd = jd + 1
        jb = jb + 1
    ja = ja + 1
    jb = 1
    jd = 1

with open('results/instab_params_matrix_3params.csv', 'w') as outfile:
    outfile.write('# Array shape: {0}\n'.format(Z.shape))
    for data_slice in Z:
        np.savetxt(outfile, data_slice, delimiter=",", fmt='%d', footer='# New slice\n')
