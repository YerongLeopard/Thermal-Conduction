import scipy as sp
import numpy as np
import matplotlib;
try:
   from PyQt4 import QtCore, QtGui
except:
   matplotlib.use("TkAgg")
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import mpl_toolkits.mplot3d.axes3d as p3

datak=pl.loadtxt('./analysis/divergencefull/6by6.dat')
dataeff=pl.loadtxt('./analysis/kappa_eff/k_eff.dat')
xx1      =[2*sp.pi/line[0] for line in datak if line[0]> 100 and line[0]<2801]
kappaSet1=[line[1] for line in datak if line[0]> 100 and line[0]<2801]
xx2      =[2*sp.pi/line[0] for line in dataeff if line[0]> 100 and line[0]<2801]
kappaSet2=[line[1] for line in dataeff if line[0]> 100 and line[0]<2801]
plt.figure(1)
plt.title("$\mathbf{Argon}$: Comparing $\kappa(q)$ and $k_\mathrm{eff}$ on Sample with Length $L$")
plt.plot(xx1,kappaSet1, '-b', label="$\kappa(\\frac{2\pi}{L})/\kappa_0$")
plt.plot(xx2,kappaSet2, '-r', label="$\kappa_\\mathrm{eff}/\kappa_0$")
pl.xlabel('$qa = 2\pi/L$', fontsize=20)
pl.ylabel('$\kappa(q)$, $\kappa_\mathrm{eff}/\kappa_0$', fontsize=20)
plt.legend(loc='best', fontsize=20)
plt.savefig("keff.png")

plt.figure(2)
plt.title("$\mathbf{Argon}$: Comparing $\kappa(q)$ and $k_\mathrm{eff}$ on Sample with Length $L$")
resSet1 = [1/element for element in kappaSet1]
resSet2 = [1/element for element in kappaSet2]
plt.plot(xx1,resSet1, '-b', label="$\kappa_0/\kappa(\\frac{2\pi}{L})$")
plt.plot(xx2,resSet2, '-r', label="$\kappa_0/\kappa_\\mathrm{eff}$")
plt.legend(loc='best', fontsize=20)
pl.xlabel('$qa = 2\pi/L$', fontsize=20)
pl.ylabel('Resistance $\kappa_0/\kappa$', fontsize=20)
plt.savefig("resistance.png")
plt.show()
# print pointSet1 # DEBUG
# print data1[0,0]# DEBUG

print('finished \n')
