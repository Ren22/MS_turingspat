import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

n = -2.
a = 0.2
d = -0.8
Du = 1.
Dv = 20.
Dm = 1000.
q = 0.001
km = 0.01
lamb = np.sqrt(Dm / km)
L = np.arange(40., 140., 1.)
l = []
M = []
kon = 0.05
koff = 0.05
laterlal_d = 20.
Mtot = 4000.
for i in L:
    # print i
    M = kon * Mtot / (laterlal_d * (kon + koff) * i)
    print(M)
    # l.append(2. * np.pi * ((i / (q*lamb**2.)) ** (n / 2.)) * (a / (2. * Du) + d / (2. * Dv)) ** (-1 / 2.)) ##control a =0.2 Du =1.0

    # l.append(2. * np.pi * M**(-n / 2.) * ((a/(2. * Du) + d / (2. * Dv)) ** (-1 / 2.)))  ##control a =0.2 Du =1.0
    l.append(2. * np.pi * ((a * M ** n) / (2. * Du) + d / (2. * Dv)) ** (-1 / 2.))  ##control a =0.2 Du =1.0

raw_data = {'size': L,
            'values': l}
df = pd.DataFrame(raw_data, columns=['size', 'values'])
df.to_csv('results/linear_resuls.csv')
# (q*lamb**2)
# *((L/(q*lamb**2))**(n/2))
# plt.title('Pattern wavelength vs system size n = 2')
plt.xlabel('L')
plt.ylabel('l')
plt.plot(L, l, '-r')
plt.show()
