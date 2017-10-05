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

# Paramters range definition
a_range = np.arange(.1, 1., .05)  # rows
b_range = np.arange(-.1, -1., -.05)  # columns
sys_size = 100.

# Matrix definition
Z = np.empty((len(a_range) + 1, len(b_range) + 1))
Z[0, 1:] = [np.round(elem, 2) for elem in b_range]  # all b_range assigned to 0 row
Z[1:, 0] = [np.round(elem, 2) for elem in a_range]  # all a_range assigned to 0 column

k = 1;
m = 1;

for i in a_range:
    for j in b_range:
        # print(i, j)
        if (i == 0): i = 1e-10
        if (j == 0): j = 1e-10
        res = criter0.subs([(a, i), (b, j), (c, 0.6), (d, -0.8),
                            (au, 5.), (u0, 0.5), (v0, 0.5), (n, -2.), (M, 100. * (1. / sys_size)),
                            (u, 0.5), (v, 0.5)])
        if res.is_Mul:
            res = 2 * Dv /Du
        Z[k, m] = res
        m = m + 1
        # print(k, m)
    k = k + 1
    m = 1

np.savetxt('results/instab_params_matrix.csv', Z, delimiter=",")

X, Y = np.meshgrid(a_range, b_range)

ax = plt.subplot(projection='3d')
ax.plot_surface(X, Y, np.transpose(Z[1:, 1:]), color='brown')  # rstride=20, cstride=20)
ax.plot_wireframe(X, Y, Dv / Du)  # Criterion surface
ax.set_xlabel('a')
ax.set_ylabel('b')
ax.set_zlabel('Crit.Value')
ax.set_title('a, b coefficienets, sys.size = {0}'.format(sys_size))

plt.show()
