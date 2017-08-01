from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def f(u, uc, kon, koff, V, L, d):
    return kon * V / (L * d) * uc - koff * u


def g(u, uc, kon, koff, V, L, d):
    return koff * L * d / V * u - kon * uc


def pde(y, t, Du, Duc, kon, koff, V, L, d, dx):
    u = y[::2]
    uc = y[1::2]

    dydt = np.empty_like(y)

    dudt = dydt[::2]
    ducdt = dydt[1::2]

    dudt[0] = f(u[0], uc[0], kon, koff, V, L, d) + Du * (-2.0 * u[0] + 2.0 * u[1]) / dx ** 2
    dudt[1:-1] = f(u[1:-1], uc[1:-1], kon, koff, V, L, d) + Du * np.diff(u, 2) / dx ** 2
    dudt[-1] = f(u[-1], uc[-1], kon, koff, V, L, d) + Du * (- 2.0 * u[-1] + 2.0 * u[-2]) / dx ** 2
    ducdt[0] = g(u[0], uc[0], kon, koff, V, L, d) + Duc * (-2.0 * uc[0] + 2.0 * uc[1]) / dx ** 2
    ducdt[1:-1] = g(u[1:-1], uc[1:-1], kon, koff, V, L, d) + Duc * np.diff(uc, 2) / dx ** 2
    ducdt[-1] = g(u[-1], uc[-1], kon, koff, V, L, d) + Duc * (-2.0 * uc[-1] + 2.0 * uc[-2]) / dx ** 2

    return dydt


h = 100.
x = np.linspace(0., 100., num=h)
print(len(x))

# #step-like function
# init_cond = [0]*(len(x))
# for i in range(int(round(0.2*len(init_cond)))):
#     init_cond[i] = 1.

# Sinusoids as initial funct
init_cond = np.empty_like(x)
init_cond_u = init_cond[::2]
init_cond_uc = init_cond[1::2]
init_cond_u[0:len(init_cond_u)] = 0.5 + 0.1 * np.cos(2. * np.pi * x[::2] / x[-1])
init_cond_uc[0:len(init_cond_u)] = 1. + 0.1 * np.cos(2. * np.pi * x[1::2] / x[-1])

t = np.linspace(0, 1000, 2000)

Du = 1.;
Duc = 100.;
kon = 0.05;
koff = 0.05;
V = 10 ** -18;
L = 10 ** -6;
d = 5. * 10 ** -6
M = 100. * (1. / x[-1]);
dx = len(x) / h
sol = odeint(pde, init_cond, t, args=(Du, Duc, kon, koff, V, L, d, dx), ml=2, mu=2)

sol_u = sol[::2]
sol_uc = sol[1::2]
# print(len(sol_u))
# print(sol_u[-1][::2])
plt.plot(x[::2], sol_u[-1][::2])
plt.plot(x[1::2], sol_uc[-1][1::2])
plt.show()
