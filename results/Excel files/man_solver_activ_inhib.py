from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def f(u, v, M, a, b, c, d, au, u0, v0, n):
    return (a * (u - u0) + b * (v - v0) - au * (u - u0) ** 3) * M ** n


def g(u, v, M, a, b, c, d, au, u0, v0, n):
    return (c * (u - u0) + d * (v - v0)) * M ** n


def pde(y, t, Du, Dv, a, b, c, d, au, u0, v0, n, M, dx):
    u = y[::2]
    v = y[1::2]

    dydt = np.empty_like(y)

    dudt = dydt[::2]
    dvdt = dydt[1::2]

    dudt[0] = f(u[0], v[0], M, a, b, c, d, au, u0, v0, n) + Du * (-2.0 * u[0] + 2.0 * u[1]) / dx ** 2
    dudt[1:-1] = f(u[1:-1], v[1:-1], M, a, b, c, d, au, u0, v0, n) + Du * np.diff(u, 2) / dx ** 2
    dudt[-1] = f(u[-1], v[-1], M, a, b, c, d, au, u0, v0, n) + Du * (- 2.0 * u[-1] + 2.0 * u[-2]) / dx ** 2
    dvdt[0] = g(u[0], v[0], M, a, b, c, d, au, u0, v0, n) + Dv * (-2.0 * v[0] + 2.0 * v[1]) / dx ** 2
    dvdt[1:-1] = g(u[1:-1], v[1:-1], M, a, b, c, d, au, u0, v0, n) + Dv * np.diff(v, 2) / dx ** 2
    dvdt[-1] = g(u[-1], v[-1], M, a, b, c, d, au, u0, v0, n) + Dv * (-2.0 * v[-1] + 2.0 * v[-2]) / dx ** 2

    return dydt


h = 100.
x = np.linspace(0., 50., num=h)
print(len(x))

# init_cond = 0.5 * np.sin(2.*np.pi*x/5.)
# plt.plot(x, init_cond)
# plt.show()

init_cond = [0] * (len(x))
for i in range(int(round(0.2 * len(init_cond)))):
    init_cond[i] = 1.
print(init_cond)

t = np.linspace(0, 1000, 2000)

Du = 1.;
Dv = 20.;
a = 0.2;
b = -0.4;
c = 0.6;
d = -0.8;
au = 5.;
u0 = 0.5;
v0 = 0.5;
n = -2;
M = 100. * (1. / x[-1]);
dx = len(x) / h
sol = odeint(pde, init_cond, t, args=(Du, Dv, a, b, c, d, au, u0, v0, n, M, dx), ml=2, mu=2)

plt.plot(x, sol[-1])
plt.show()
