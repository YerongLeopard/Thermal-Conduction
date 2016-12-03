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
data=pl.loadtxt('ABC.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[1])
  b.append(line[2])
  c.append(line[3])
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'ob',markersize=8)
pl.plot(x1,a,'-b',markersize=8)
pl.xlabel('$1/N_x^2$',fontsize=25)
pl.ylabel('$A$',fontsize=25)
pl.figure(2)
pl.plot(x2,b,'ob',markersize=8)
pl.plot(x2,b,'-b',markersize=8)
pl.xlabel('$1/N_x$',fontsize=25)
pl.ylabel('$B$',fontsize=25)
pl.figure(3)
pl.plot(x3,c,'ob',markersize=8)
pl.plot(x3,c,'-b',markersize=8)
pl.xlabel('$N_x^2$',fontsize=25)
pl.ylabel('$C$',fontsize=25)

#########################################50
data=pl.loadtxt('ABC50.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[1])
  b.append(line[2])
  c.append(line[3])
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'oy',markersize=8,label="N_z>50")
pl.plot(x1,a,'-y',markersize=8)

pl.figure(2)
pl.plot(x2,b,'oy',markersize=8,label="N_z>50")
pl.plot(x2,b,'-y',markersize=8)
pl.figure(3)
pl.plot(x3,c,'oy',markersize=8,label="N_z>50")
pl.plot(x3,c,'-y',markersize=8)

'''
########################################400
data=pl.loadtxt('ABC400.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[1])
  b.append(line[2])
  c.append(line[3])
  x1.append(1.0/line[0]**2)
  x2.append(1.0/line[0])
  x3.append(line[0]**2)
pl.figure(1)
pl.plot(x1,a,'og',markersize=8,label="N_z>400")
pl.plot(x1,a,'-g',markersize=8)

pl.figure(2)
pl.plot(x2,b,'og',markersize=8,label="N_z>400")
pl.plot(x2,b,'-g',markersize=8)
pl.figure(3)
pl.plot(x3,c,'og',markersize=8,label="N_z>400")
pl.plot(x3,c,'-g',markersize=8)


########################################800
data=pl.loadtxt('ABC800.dat')
a=[]
b=[]
c=[]
x1=[]
x2=[]
x3=[]
for line in data:
  a.append(line[1])
  b.append(line[2])
  c.append(line[3])
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
  a.append(line[1])
  b.append(line[2])
  c.append(line[3])
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

pl.legend(loc='best')


pl.show()
