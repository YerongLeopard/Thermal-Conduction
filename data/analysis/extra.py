import matplotlib;
try:
   from PyQt4 import QtCore, QtGui
except:
   matplotlib.use("TkAgg")
import scipy as sp
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import mpl_toolkits.mplot3d.axes3d as p3
import math as m
data=pl.loadtxt('ABC10.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[[1]]*3)
  b.append(line[[2]]*3)
  c.append(line[[3]]*3)
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'ob',markersize=8,label="$N_x>10$")
pl.plot(x1,a,'-b',markersize=8)
pl.xlabel('$1/N_y^2$',fontsize=25)
h=pl.ylabel('$A$',fontsize=25)
h.set_rotation(0)
pl.legend(loc='best')

pl.figure(2)
pl.plot(x2,b,'ob',markersize=8,label="$N_x>10$")
pl.plot(x2,b,'-b',markersize=8)
pl.xlabel('$1/N_y$',fontsize=25)
h=pl.ylabel('$B$',fontsize=25)
h.set_rotation(0)
pl.legend(loc='best')

pl.figure(3)
pl.plot(x3,c,'ob',markersize=8,label="$N_x>10$")
pl.plot(x3,c,'-b',markersize=8)
pl.xlabel('$N_y^2$',fontsize=25)
h=pl.ylabel('$C$',fontsize=25)
h.set_rotation(0)
pl.legend(loc='best')
#########################################50
data=pl.loadtxt('ABC50.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[[1]]*3)
  b.append(line[[2]]*3)
  c.append(line[[3]]*3)
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'oy',markersize=8,label="$N_x>50$")
pl.plot(x1,a,'-y',markersize=8)
pl.legend(loc='best')

ax = plt.gca()
# Customize tick marks
ax.tick_params(which='both',width=2,length=8) # setting length and width of tick marks

# Customize tick labels
for idx, tick in enumerate(ax.xaxis.get_major_ticks()):
    tick.label1.set_fontsize(12)
    tick.label1.set_fontweight('bold')

for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(12)
    tick.label1.set_fontweight('bold')

pl.legend(loc='best')
pl.figure(2)
pl.plot(x2,b,'oy',markersize=8,label="$N_x>50$")
pl.plot(x2,b,'-y',markersize=8)
pl.legend(loc='best')

ax = plt.gca()
# Customize tick marks
ax.tick_params(which='both',width=2,length=8) # setting length and width of tick marks

# Customize tick labels
for idx, tick in enumerate(ax.xaxis.get_major_ticks()):
    tick.label1.set_fontsize(12)
    tick.label1.set_fontweight('bold')

for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(12)
    tick.label1.set_fontweight('bold')

pl.legend(loc='best')


pl.figure(3)
pl.plot(x3,c,'oy',markersize=8,label="$N_x>50$")
pl.plot(x3,c,'-y',markersize=8)
pl.legend(loc='best')

ax = plt.gca()
# Customize tick marks
ax.tick_params(which='both',width=2,length=8) # setting length and width of tick marks

# Customize tick labels
for idx, tick in enumerate(ax.xaxis.get_major_ticks()):
    tick.label1.set_fontsize(12)
    tick.label1.set_fontweight('bold')

for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(12)
    tick.label1.set_fontweight('bold')

pl.legend(loc='best')
########################################400
data=pl.loadtxt('ABC400.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[[1]]*3)
  b.append(line[[2]]*3)
  c.append(line[[3]]*3)
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'or',markersize=8,label="$N_x>400$")
pl.plot(x1,a,'-r',markersize=8)
pl.legend(loc='best')

pl.figure(2)
pl.plot(x2,b,'or',markersize=8,label="$N_x>400$")
pl.plot(x2,b,'-r',markersize=8)
pl.legend(loc='best')

pl.figure(3)
pl.plot(x3,c,'or',markersize=8,label="$N_x>400$")
pl.plot(x3,c,'-r',markersize=8)
pl.legend(loc='best')

'''
########################################800
data=pl.loadtxt('ABC800.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[[1,2,3]]*3)
  b.append(line[[1,2,3]]*3)
  c.append(line[[1,2,3]]*3)
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'or',markersize=8,label="N_z>800")
pl.plot(x1,a,'-r',markersize=8)

pl.figure(2)
pl.plot(x2,b,'or',markersize=8,label="N_z>800")
pl.plot(x2,b,'-r',markersize=8)
pl.figure(3)
pl.plot(x3,c,'or',markersize=8,label="N_z>800")
pl.plot(x3,c,'-r',markersize=8)

########################################1000
data=pl.loadtxt('ABC1000.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[[1,2,3]]*3)
  b.append(line[[1,2,3]]*3)
  c.append(line[[1,2,3]]*3)
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'oc',markersize=8,label="N_z>1000")
pl.plot(x1,a,'-c',markersize=8)
pl.legend(loc='best')
pl.figure(2)
pl.plot(x2,b,'oc',markersize=8,label="N_z>1000")
pl.plot(x2,b,'-c',markersize=8)
pl.legend(loc='best')
pl.figure(3)
pl.plot(x3,c,'oc',markersize=8,label="N_z>1000")
pl.plot(x3,c,'-c',markersize=8)
'''


pl.show()
