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
kappaSet1=[line[1] for line in data1 if line[0]> 99       and line[0]<2001         and 0==line[0]%10       ]
kappaSet2=[line[1] for line in data2 if int(line[0]/4)>99 and int(line[0]/4)<2001  and 0==int(line[0]/4)%10]
xx = [2*sp.pi/length for length in np.arange(100,2001,10)]
plt.plot(xx,kappaSet1, '--b', label="$\kappa(\\frac{2\pi}{L})/\kappa_0$ on sample with length $L$")
plt.plot(xx,kappaSet2,  '.r', label="$\kappa(\\frac{2\pi}{L})/\kappa_0$ on sample with length $4L$")
pl.xlabel('$qa = 2\pi/L$', fontsize=20)
pl.ylabel('$\kappa(q=\\frac{2\pi}{L})/\kappa_0$', fontsize=20)
plt.legend(loc='best')
plt.savefig('comparing_q_with_4q.png')
plt.show()
# print pointSet1 # DEBUG
# print data1[0,0]# DEBUG
print('finished \n')
