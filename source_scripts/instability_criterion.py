import numpy as np

L = 60.
d = 20.
V = 2700.
kon = 0.05
koff = 0.05
Du = 100.
Dv = 1.
d = Dv / Du

fu = -koff
fv = kon * V / (L * d)
gu = koff * L * d / V
gv = -kon

# print('fu*gv={0}, fv*gu = {1}, fu*gv-fv*gu = {2}'.format(fu*gv,fv*gu, fu*gv-fv*gu))
# print('np.sqrt(fu*gv-fv*gu) = {0}'.format(np.sqrt(fu*gv-fv*gu)))
# print('np.sqrt(-fv*gu)= {0}'.format(np.sqrt(-fv*gu)))
# print('np.sqrt(fu*gv-fv*gu)+np.sqrt(-fv*gu)={0}'.format(np.sqrt(fu*gv-fv*gu)+np.sqrt(-fv*gu)))
#
# print('final_res={0}'.format((np.sqrt(round(fu*gv-fv*gu))+np.sqrt(-fv*gu))**2/fu**2))

print('1) fu+gv < 0 = {0}, '
      '2) fu*gv-fv*gu > 0 = {1}, '
      '3) dfu+gv > 0 = {2}, '
      '4)(dfu+gv)^2-4d(fugv-fvgu) > 0 = {3}'.
      format(fu + gv, fu * gv - fv * gu, d * fu + gv, (d * fu + gv) ** 2 - 4 * d * (fu * gv - fv * gu)))
