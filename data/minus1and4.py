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

data1=pl.loadtxt('./analysis/divergencefull/6by6.dat')
data2=pl.loadtxt('./analysis/divergence-4/6by6-4.dat')
kappaSet1=[line[1] for line in data1 if line[0]>399 and line[0]<12001 and 0==line[0]%4]
kappaSet4=[line[1] for line in data2 if line[0]>399 and line[0]<12001 and 0==line[0]%4]
xx = [2*sp.pi/length for length in np.arange(400,12001,4)]
yy = [2*kappa_4 - kappa_1 for kappa_1, kappa_4 in zip(kappaSet1, kappaSet4)]
plt.plot(xx, yy, '-b', label="$(2*\kappa(4q)-\kappa(q))/\kappa_0$ on sample with length $L$")
pl.xlabel('$qa = 2\pi/L$', fontsize=20)
pl.ylabel('$(2*\kappa(4q)-\kappa(q))/\kappa_0$', fontsize=20)
plt.legend(loc='best')
plt.savefig('minus.png')
plt.show()
# print pointSet1 # DEBUG
# print data1[0,0]# DEBUG
print('finished \n')
