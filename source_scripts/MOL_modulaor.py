from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

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
init_cond_uc[0:len(init_cond_u)] = 0.5 + 0.1 * np.cos(2. * np.pi * x[1::2] / x[-1])

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

##3D plot (activator)
ax_u = plt.subplot(223, projection='3d')
ax_uc = plt.subplot(224, projection='3d')

# Grab data for u, uc
xx_u, xx_uc = x[::2], x[1::2]
yy_u = t;
yy_uc = t
zz_u, zz_uc = [], []
for i in sol:
    zz_u.append(i[::2])
    zz_uc.append(i[1::2])
XX_u, YY_u = np.meshgrid(xx_u, yy_u);
XX_uc, YY_uc = np.meshgrid(xx_uc, yy_uc)
ZZ_u = zz_u;
ZZ_uc = zz_uc

# Plot a basic wireframe for u
ax_u.plot_surface(XX_u, YY_u, ZZ_u, rstride=20, cstride=20)
ax_u.set_xlabel('Space')
ax_u.set_ylabel('Time')
ax_u.set_zlabel('Value')
ax_u.set_title('Activator profile')

# Plot a basic wireframe for uc
ax_uc.plot_surface(XX_uc, YY_uc, ZZ_uc, rstride=20, cstride=20)
ax_uc.set_xlabel('Space')
ax_uc.set_ylabel('Time')
ax_uc.set_zlabel('Value')
ax_uc.set_title('Activator profile')
# Show plot

plt.show()