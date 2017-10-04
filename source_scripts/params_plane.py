import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from sympy import *
import sympy

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
a_range = np.arange(0.1, 0.3, .1)
b_range = np.arange(-0.3, -1., -.1)

sys_size = 50.
Z = np.empty((len(a_range), len(b_range)));
k = 0;
m = 0;

for i in a_range:
    for j in b_range:
        # print(i, j)
        if (i == 0): i = 1e-10
        if (j == 0): j = 1e-10
        res = criter0.subs([(a, i), (b, j), (c, 0.6), (d, -0.8),
                            (au, 5.), (u0, 0.5), (v0, 0.5), (n, -2.), (M, 100. * (1. / sys_size)),
                            (u, 0.5), (v, 0.5)])
        if res.is_Mul:
            res = 0.
        Z[k][m] = res
        m = m + 1
        # print(k, m)
    k = k + 1
    m = 0

X, Y = np.meshgrid(a_range, b_range)
print(Z)

ax = plt.subplot(projection='3d')
ax.plot_surface(X, Y, Z)  # rstride=20, cstride=20)
ax.plot_surface(X, Y, Dv / Du)  # Criterion surface
ax.set_xlabel('a')
ax.set_ylabel('b')
ax.set_zlabel('Crit.Value')
ax.set_title('a, b coefficienets, sys.size = {0}'.format(sys_size))

plt.show()
