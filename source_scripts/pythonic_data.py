from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

xl = pd.ExcelFile('../results/Excel files/Compare_selected_data.xlsx')
df = xl.parse('Sheet1')
columns_nums = len(df.columns)
df = df.values

size_list = []
wv_list = []

for i in np.arange(1, columns_nums):
    size_list = np.append(size_list, df[:, i])

ss = [df[:, 0]] * (columns_nums - 1)

for j in np.arange(0, len(ss)):
    wv_list = np.append(wv_list, ss[j])

font = {'family': 'normal',
        'weight': 'bold',
        'size': 14}

matplotlib.rc('font', **font)
plt.hist2d(wv_list, size_list, bins=20, cmap='Greens')
plt.xlabel('System size')
plt.ylabel('Pattern wavelength')
plt.colorbar()

plt.show()
