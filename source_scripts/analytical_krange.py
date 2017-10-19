import numpy as np
from sympy import *
import matplotlib.pyplot as plt

a, b, c, d, au, u0, v0, n, M, u, v = symbols('a, b, c, d, au, u0, v0, n, M, u, v')

f = a * (u - u0) * M ** n + b * (v - v0) - (au * (u - u0) ** 3)
g = c * (u - u0) + d * (v - v0)
fu = diff(f, u)
fv = diff(f, v)
gu = diff(g, u)
gv = diff(g, v)

Du = 1.
Dv = 90.

k_range = np.arange(0., 100., .1)
# gamma = size**2/Du
for size in np.arange(5., 120., 1.):
    h_arr = []
    params_list = [(a, 0.4), (b, -2.), (c, 2.), (d, -2.),
                   (au, 5.), (u0, 0.5), (v0, 0.5), (n, -2.),
                   (u, 0.5), (v, 0.5), (M, 20. * (1. / size))]
    for k in k_range:
        # params_list.append((M, 20. * (1. / size)))
        h = (Dv / Du) * k ** 4 - (size ** 2 / Du) * ((Dv / Du) * fu + gv) * k ** 2 + \
            (size ** 2 / Du) * (fu * gv - fv * gu)
        h_arr.append(h.subs(params_list))
        if len(h_arr) == len(k_range):
            # h_arr_pos = h_arr
            # h_arr_neg = h_arr
            # h_arr_pos[h_arr_pos <= 0] = np.nan
            # h_arr_neg[h_arr_neg > 0] = np.nan
            plt.xlabel('k')
            plt.ylabel('h')
            plt.title('h=h(k), size={0}'.format(size))
            plt.plot(k_range ** 2, h_arr, 'g')
            plt.savefig('results/h_vs_k/size_{0}.png'.format(size))
            plt.close()
            # plt.show()
