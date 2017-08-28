import numpy as np
import matplotlib.pyplot as plt

h = 60.
x = np.linspace(0., 60., num=h)

# init_cond = np.empty_like(x)
# # init_cond_u = init_cond[::2]
# # init_cond_v = init_cond[1::2]
# u_02 = int(round(0.2 * len(init_cond)))
# u_08 = int(round(0.8 * len(init_cond)))
# init_cond[0:u_02] = [0.6]*u_02
# init_cond[u_02:u_08] = [0.5]*(u_08-u_02)
# init_cond[u_08:len(init_cond)] = [0.6] * (len(init_cond)-u_08)

# init_cond = np.empty_like(x)
# init_cond_u = init_cond[::2]
# init_cond_v = init_cond[1::2]
# init_cond_u[0:len(init_cond_u)] = 0.5 + 0.1 * np.cos(2. * np.pi * x[::2] / x[-1])
# init_cond_v[0:len(init_cond_u)] = 1. + 0.1 * np.cos(2. * np.pi * x[1::2] / x[-1])

##FFT of init cond

u_fft_y = init_cond - np.mean(init_cond)
u_fft_x = x
u_fft_x_norm = x[-1] * np.array(u_fft_x)

Y1 = np.fft.fft(u_fft_y)
N = len(Y1) / 2 + 1
dt = x[-1] / h
print (dt)
fa = 1.0 / dt
X = np.linspace(0, fa / 2, N, endpoint=True)

plt.plot(X, 2.0 * np.abs(Y1[:N]) / N)
plt.show()
