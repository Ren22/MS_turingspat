from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

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
x = np.linspace(0., 19., num=h)

##Sinusoids as initial function
# init_cond = np.empty_like(x)
# init_cond_u = init_cond[::2]
# init_cond_v = init_cond[1::2]
# init_cond_u[0:len(init_cond_u)] = 0.5 + 0.1 * np.cos(2. * np.pi * x[::2] / x[-1])
# init_cond_v[0:len(init_cond_u)] = 1. + 0.1 * np.cos(2. * np.pi * x[1::2] / x[-1])

##step like as initial_function
init_cond = np.empty_like(x)
init_cond_u = init_cond[::2]
init_cond_v = init_cond[1::2]
u_02 = int(round(0.2 * len(init_cond_u)));
v_02 = int(round(0.2 * len(init_cond_v)))
u_08 = int(round(0.8 * len(init_cond_u)));
v_08 = int(round(0.8 * len(init_cond_v)))
init_cond_u[0:u_02] = [0.6] * u_02
init_cond_u[u_02:u_08] = [0.5] * (u_08 - u_02)
init_cond_u[u_08:len(init_cond_u)] = [0.6] * (len(init_cond_u) - u_08)
init_cond_v[0:v_02] = [0.6] * v_02
init_cond_v[u_02:u_08] = [0.5] * (v_08 - v_02)
init_cond_v[v_08:len(init_cond_v)] = [0.6] * (len(init_cond_v) - v_08)

t = np.linspace(0, 5000, 5000)

Du = 1.;
Dv = 20.;
a = 0.2;
b = -0.4;
c = 0.6;
d = -0.8;
au = 5.;
u0 = 0.5;
v0 = 0.5;
n = -2.;
M = 100. * (1. / x[-1]);
dx = len(x) / h
sol = odeint(pde, init_cond, t, args=(Du, Dv, a, b, c, d, au, u0, v0, n, M, dx), ml=2, mu=2)

##2D plot
plt.figure(figsize=(10, 7))
# activator
plt.subplot(221)
plt.title('Activator profile')
plt.xlabel('System size')
plt.ylabel('Concentration')
plt.plot(x[::2], sol[-1][::2], '-b')

# inhibitor
plt.subplot(222)
plt.title('Inhibitor profile')
plt.xlabel('System size')
plt.ylabel('Concentration')
plt.plot(x[1::2], sol[-1][1::2], '-r')

## Modulator
plt.subplot(223)
plt.title('Modulator profile')
plt.xlabel('System size')
plt.ylabel('Concentration')
plt.plot(x, [M] * len(x), '-g')

##3D plot (activator)
ax = plt.subplot(224, projection='3d')

# Grab data.
xx = x[::2]
yy = t
zz = []
for i in sol:
    zz.append(i[::2])

XX, YY = np.meshgrid(xx, yy);
ZZ = zz

# Plot a basic wireframe
ax.plot_surface(XX, YY, ZZ, rstride=20, cstride=20)
ax.set_xlabel('Space')
ax.set_ylabel('Time')
ax.set_zlabel('Value')
ax.set_title('Activator profile')

# Show plot

plt.show()