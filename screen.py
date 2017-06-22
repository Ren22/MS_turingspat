import os
import numpy as np

#au = 0.2 # activator kinetic constant
#Du = 1 # activator diffusion constant
time_period = 50
times = 30

for au in np.arange(0.17, 0.26, 0.01):
    for Du in np.arange(0.7, 1.3, 0.1):
        for size in np.arange(40, 121, 5):
            size_segments = 2 * size
            os.system("python wavesearch.py {0} {1} {2} {3} {4} {5}".format(size, size_segments, time_period, times, au, Du))