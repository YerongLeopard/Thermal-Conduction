import scipy as sp
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import mpl_toolkits.mplot3d.axes3d as p3
import math as m
data=pl.loadtxt('cross1.dat')
x=[]
y=[]
for line in data:
  x.append(m.log(line[0],10))
  y.append(line[1])
print y
pl.plot(x,y,'.',markersize=10)
pl.xlabel('Width of Cross Section $N_1$',fontsize=25)
pl.ylabel('$\kappa_{N_1,N_3}$',fontsize=25)
plt.xticks([i for i in range(0,2+1,1)],[10**i for i in range(0,2+1,1)],fontsize=14)
pl.yticks(fontsize=14)
pl.show()
