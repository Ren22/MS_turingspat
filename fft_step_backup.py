import numpy as np
import matplotlib.pyplot as plt

# Grid
m = 10.
segm = m * 5.
x = np.linspace(0, m, segm, endpoint=False)  # np.arange(0, m)
k = .001  # diff const

# Time steps
tmax = 15.
dt = 0.05

# Initial data
u = 5 * np.sin(x / 2 * np.pi)
u_fft = np.fft.fft(u)
N = len(u_fft)
dx = m / segm
fa = 1 / dx  # sample frequency
xi = np.linspace(0, fa, N, endpoint=False)

for time in np.arange(0., tmax, dt):
    E = np.exp(-k * time * xi ** 2.)
    u_fft = E * u_fft

    plt.clf()
    u = np.squeeze(np.real(np.fft.ifft(u_fft)))
    plt.plot(x, u)
    plt.xlim((0, m))
    plt.ylim((-10., 10.))
    plt.title('t=' + str(time))
    plt.show()
