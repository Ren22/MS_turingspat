import numpy as np
import matplotlib.pyplot as plt

a = 0.2
b = -0.4
c = 0.6
d = -0.8
au = 5.
u0 = 0.5
v0 = 0.5

u = np.arange(-.2, 1, .01)
v1 = []
v2 = []
for i in u:
    v1.append((1 / b) * (au * (i - u0) ** 3 - a * (i - u0)) + v0)
    v2.append((-c / d) * (i - u0) + v0)
plt.title('Phase portrait')
plt.xlabel('u')
plt.ylabel('v')
plt.plot(u, v1, '-r')
plt.plot(u, v2, '-b')
plt.show()
