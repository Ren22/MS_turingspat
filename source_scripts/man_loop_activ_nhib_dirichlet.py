from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import axes3d

Du = 1.
Dv = 20.
a = 0.2
b = -0.4
c = 0.6
d = -0.8
au = 5.
u0 = 0.5
v0 = 0.5
n = -2.

def f(u, v, M, a, b, c, d, au, u0, v0, n):
    return a * (u - u0) * M ** n + b * (v - v0) - (au * (u - u0) ** 3)

def g(u, v, M, a, b, c, d, au, u0, v0, n):
    return c * (u - u0) * M ** n + d * (v - v0)

def pde(y, t, Du, Dv, a, b, c, d, au, u0, v0, n, M, dx):
    u = y[::2]
    v = y[1::2]

    dydt = np.empty_like(y)
    dudt = dydt[::2]
    dvdt = dydt[1::2]

    # # Dirichlet
    # dudt[0] = 0
    # dudt[1:-1] = f(u[1:-1], v[1:-1], M, a, b, c, d, au, u0, v0, n) + Du * np.diff(u, 2) / dx ** 2
    # dudt[-1] = 0
    # dvdt[0] = 0
    # dvdt[1:-1] = g(u[1:-1], v[1:-1], M, a, b, c, d, au, u0, v0, n) + Dv * np.diff(v, 2) / dx ** 2
    # dvdt[-1] = 0

    # # Neumann
    dudt[0] = f(u[0], v[0], M, a, b, c, d, au, u0, v0, n) + Du * (-2.0 * u[0] + 2.0 * u[1]) / dx ** 2
    dudt[1:-1] = f(u[1:-1], v[1:-1], M, a, b, c, d, au, u0, v0, n) + Du * np.diff(u, 2) / dx ** 2
    dudt[-1] = f(u[-1], v[-1], M, a, b, c, d, au, u0, v0, n) + Du * (- 2.0 * u[-1] + 2.0 * u[-2]) / dx ** 2
    dvdt[0] = g(u[0], v[0], M, a, b, c, d, au, u0, v0, n) + Dv * (-2.0 * v[0] + 2.0 * v[1]) / dx ** 2
    dvdt[1:-1] = g(u[1:-1], v[1:-1], M, a, b, c, d, au, u0, v0, n) + Dv * np.diff(v, 2) / dx ** 2
    dvdt[-1] = g(u[-1], v[-1], M, a, b, c, d, au, u0, v0, n) + Dv * (-2.0 * v[-1] + 2.0 * v[-2]) / dx ** 2

    return dydt

# Main Solver
t = np.linspace(0, 5000, 5000)

for size in np.arange(10., 120., 1.):
    discretizing_factor = 6.
    size_segments = discretizing_factor * size
    x = np.linspace(0., size, size_segments)
    M = 100. * (1. / x[-1])
    dx = (size / size_segments) * 2.


    ## Cosinusoids as initial function
    # init_cond = np.empty_like(x)
    # init_cond_u = init_cond[::2]
    # init_cond_v = init_cond[1::2]
    # init_cond_u[0:len(init_cond_u)] = 0.5 + 0.1 * np.cos(2. * np.pi * x[::2] / x[-1])
    # init_cond_v[0:len(init_cond_u)] = 1. + 0
    # .1 * np.cos(2. * np.pi * x[1::2] / x[-1])

    ## Step like as initial_function
    # init_cond = np.empty_like(x)
    # init_cond_u = init_cond[::2]
    # init_cond_v = init_cond[1::2]
    # u_045 = int(round(0.45 * len(init_cond_u)))
    # v_045 = int(round(0.45 * len(init_cond_v)))
    # u_055 = int(round(0.55 * len(init_cond_u)))
    # v_055 = int(round(0.55 * len(init_cond_v)))
    # init_cond_u[0:u_045] = [0.1] * u_045
    # init_cond_u[u_045:u_055] = [0.15] * (u_055 - u_045)
    # init_cond_u[u_055:len(init_cond_u)] = [0.1] * (len(init_cond_u) - u_055)
    # init_cond_v[0:v_045] = [0.1] * v_045
    # init_cond_v[u_045:u_055] = [0.1] * (v_055 - v_045)
    # init_cond_v[v_055:len(init_cond_v)] = [0.1] * (len(init_cond_v) - v_055)

    # Step like as initial_function
    init_cond = np.empty_like(x)
    init_cond_u = init_cond[::2]
    init_cond_v = init_cond[1::2]
    u_02 = int(round(0.2 * len(init_cond_u)))
    v_02 = int(round(0.2 * len(init_cond_v)))
    u_08 = int(round(0.8 * len(init_cond_u)))
    v_08 = int(round(0.8 * len(init_cond_v)))
    init_cond_u[0:u_02] = [0.52] * u_02
    init_cond_u[u_02:u_08] = [0.4] * (u_08 - u_02)
    init_cond_u[u_08:len(init_cond_u)] = [0.52] * (len(init_cond_u) - u_08)
    init_cond_v[0:v_02] = [0.52] * v_02
    init_cond_v[u_02:u_08] = [0.4] * (v_08 - v_02)
    init_cond_v[v_08:len(init_cond_v)] = [0.52] * (len(init_cond_v) - v_08)

    sol = odeint(pde, init_cond, t, args=(Du, Dv, a, b, c, d, au, u0, v0, n, M, dx), ml=2, mu=2)

    # Activator FFT source files
    u_fft_y = sol[-1][::2] - np.mean(sol[-1][::2])
    u_fft_x = x[::2]

    # Activator FFT calculation
    Y1 = np.fft.fft(u_fft_y)
    N = len(Y1) / 2 + 1
    # dt = x[-1] / size_segments
    fa = 1.0 / (dx)
    X = np.linspace(0, fa / 2, N, endpoint=True)

    # fig = plt.figure
    plt.plot(X, 2.0 * np.abs(Y1[:N] / N))
    plt.savefig('results/plots/Fourier_{0}_{1}_{2}.png'.format(size, a, Du))
    plt.close()

    if (len(X) == len(2.0 * np.abs(Y1[:N] / N))):
        u_maxx = (np.argmax(2.0 * np.abs(Y1[:N] / N)))  # rename u_maxx
        wavelen = np.around(1. / X[u_maxx])
    # print(wavelen * 2.)

    # Write out output to the file
    df_new = pd.DataFrame([[size, size_segments, a, Du, wavelen]],
                          columns=['size', 'size_segments', 'a', 'Du', 'wavelength'])

    try:
        df = pd.read_csv('results/data.csv')
        df = df.append(df_new, ignore_index=True)
    except:
        df = df_new

    df.to_csv('results/data.csv', index=False)

    # 2D plot
    fig = plt.figure(figsize=(11.5, 7.5))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2, 2])

    # Activator
    ax0 = plt.subplot(gs[0])
    ax0.set_title('Activator profile')
    ax0.set_xlabel('System size')
    ax0.set_ylabel('Concentration')
    ax0.plot(x[::2], sol[-1][::2], '-b')

    # Inhibitor
    ax1 = plt.subplot(gs[1])
    ax1.set_title('Inhibitor profile')
    ax1.set_xlabel('System size')
    ax1.set_ylabel('Concentration')
    ax1.plot(x[1::2], sol[-1][1::2], '-r')

    # Modulator
    ax2 = plt.subplot(gs[2])
    ax2.set_title('Modulator profile')
    ax2.set_xlabel('System size')
    ax2.set_ylabel('Concentration')
    ax2.plot(x, [M] * len(x), '-g')
    plt.tight_layout()

    # 3D plot (activator)
    ax = plt.subplot(224, projection='3d')

    # Grab data for 3D
    xx = u_fft_x
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
    # plt.show()

    # Save Plot
    plt.savefig('results/plots/{0}_{1}_{2}_{3}.png'.format(size, size_segments, a, Du))

    # Close Plot

    plt.close()