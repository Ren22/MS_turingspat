from neuron import h, rxd
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d
start_time1 = time.time()

size = float(sys.argv[1])
size_segments = int(sys.argv[2])
time_period = float(sys.argv[3])
times = int(sys.argv[4])
au = float(sys.argv[5])
Du = float(sys.argv[6])

def plot_it(time_period, times, u_timespace):
    # Activator plt
    y1 = u.nodes.concentration
    x1 = u.nodes.x

    # convert x from normalized position to microns
    x1 = dend.L * np.array(x1)

    plt.figure(figsize=(20, 10))
    plt.subplot(221)
    plt.title('Activator profile')
    plt.xlabel('System size')
    plt.ylabel('Concentration')
    plt.plot(x1, y1, '-b')

    # Inhibitor plt
    y2 = z.nodes.concentration
    x2 = z.nodes.x

    # convert x from normalized position to microns
    x2 = dend.L * np.array(x2)

    plt.subplot(222)
    plt.title('Inhibitor profile')
    plt.xlabel('System size')
    plt.ylabel('Concentration')
    plt.plot(x2, y2, '-r')

    # Modulator plt
    y3 = v.nodes.concentration
    x3 = v.nodes.x

    # convert x from normalized position to microns
    x3 = dend.L * np.array(x3)

    plt.subplot(223)
    plt.title('Modulator profile')
    plt.xlabel('System size')
    plt.ylabel('Concentration')
    plt.plot(x3, y3, '-g')

    #3D plot
    ax = plt.subplot(224, projection='3d')

    # Grab data.
    xx = u_fft_x_norm
    yy = [i*time_period for i in np.arange(0, times)]
    zz = u_timespace

    XX, YY = np.meshgrid(xx, yy)
    ZZ = zz

    # Plot a basic wireframe.
    ax.plot_surface(XX, YY, ZZ, rstride=20, cstride=20)
    ax.set_xlabel('Space')
    ax.set_ylabel('Time')
    ax.set_zlabel('Value')
    ax.set_title('Activator profile')

    #plt.show() Interactive 3D mode

def wave_search(size, size_segments, time_period, times, au, Du):
    # needed for standard run system
    h.load_file('stdrun.hoc')

    global dend

    dend = h.Section()
    dend.L = size
    dend.nseg = size_segments

    # WHERE the dynamics will take place
    where = rxd.Region(h.allsec())

    # WHO the actors are
    global u
    global z
    global v

    u = rxd.Species(where, d=Du, initial=0.5)  # activator
    z = rxd.Species(where, d=20, initial=0.5)  # inhibitor
    v = rxd.Species(where, d=0, initial=(1 / dend.L) * 30)  # modulator

    # HOW they act

    a = au;
    b = -0.4;
    c = 0.6;
    d = -0.8;
    u0 = 0.5;
    z0 = 0.5;
    av = 5.0;
    kz = 0.001;

    bistable_reaction1 = rxd.Rate(u, (a * (u - u0) + b * (z - z0) - av * (u - u0) ** 3) * (v ** -2))
    bistable_reaction2 = rxd.Rate(z, (c * (u - u0) + d * (z - z0)) * (v ** -2))

    # initial conditions
    h.finitialize()
    for node in u.nodes:
        if node.x < .2: node.concentration = 0.6
        if node.x > .8: node.concentration = 0.6
    for node in z.nodes:
        if node.x < .2: node.concentration = 0.6
        if node.x > .8: node.concentration = 0.6

    # Setting up time frame
    global u_timespace

    T_d = times
    T = time_period
    u_timespace = []

    for i in np.arange(0, T_d):
        h.continuerun(i * T)
        u_timespace.append(u.nodes.concentration)

    # activator FFT source files
    u_fft_y = u.nodes.concentration
    u_fft_y = u_fft_y - np.mean(u_fft_y)
    u_fft_x = u.nodes.x
    global u_fft_x_norm
    u_fft_x_norm = dend.L * np.array(u_fft_x)

    # inhibitor FFT source files
    z_fft_y = z.nodes.concentration
    z_fft_y = z_fft_y - np.mean(z_fft_y)
    z_fft_x = z.nodes.x
    z_fft_x_norm = dend.L * np.array(u_fft_x)

    # activator FFT
    Y1 = np.fft.fft(u_fft_y)
    N = len(Y1) / 2 + 1
    dt = dend.L / dend.nseg
    fa = 1.0 / dt
    X = np.linspace(0, fa / 2, N, endpoint=True)

    # inhibitor FFT
    Y2 = np.fft.fft(z_fft_y)
    X2 = np.linspace(0, fa / 2, N, endpoint=True)
    #
    # if ((np.amax(Y1) - np.amin(Y1) < .01) or (np.amax(Y2) - np.amin(Y2) < .01)):
    #     return 0

    if (len(X) == len(2.0 * np.abs(Y1[:N] / N))):
        u_maxx = (np.argmax(2.0 * np.abs(Y1[:N] / N)))
        wavelen = np.around(1 / X[u_maxx])

    plot_it(time_period, times, u_timespace)
    plt.savefig('results/plots/{0}_{1}_{2}_{3}_{4}_{5}.png'.format(size, size_segments, time_period, times, au, Du))

    return wavelen

wl = wave_search(size, size_segments, time_period, times, au, Du)

# if (wl == 0):


runtime = time.time() - start_time1

df_new = pd.DataFrame([[size, size_segments, time_period, times, au, Du, wl, runtime]],
                      columns=['size', 'size_segments', 'time_period', 'times', 'a_u', 'Du', 'wl', 'runtime'])

try:
    df = pd.read_csv('results/data.csv')
    df = df.append(df_new, ignore_index=True)
except:
    df = df_new

df.to_csv('results/data.csv', index=False)
print(df)