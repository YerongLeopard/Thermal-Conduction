import numpy as np
import scipy.optimize as optimization
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
import scipy.optimize as optimization
import math as m
PI=3.1415926535898
data=pl.loadtxt('parameters.dat')
e=data[0]
LENX=data[1]
LENZ=data[2]
SLAB=int(data[3])
scale=LENZ/SLAB
x1=np.arange(1000*(SLAB-1))
x1=x1/1000.0
def func(x,a,b): 
 return (-a*np.cos(2*PI*(x+1)/SLAB)+b) 
x =np.arange(SLAB)
y= pl.loadtxt('slab.dat')
err=pl.loadtxt('error.dat')
y=y[-1]
err=err[-1]
j=0
x0=np.array([0.0,0.0])
para,corv=optimization.curve_fit(func,x,y,x0,err)
T=para[0]
T0=para[1]
kappa=e*scale/T/4/(np.sin(PI/SLAB)**2)/(LENX**2)
Terror=np.sqrt(corv[0][0])
kerror=Terror/T*kappa

plt.scatter(x,y,s=10)
#calculation
plt.plot(x1,func(x1,T,T0),'-r',label=str(SLAB)+'$\kappa$='+str(round(kappa,4))+'\n'+str(round(kerror,4))+'\n'+str(round(kerror/kappa,4)))
plt.xticks(np.arange(SLAB))
plt.errorbar(x,y,yerr=(err),linestyle="None")

pl.xlabel('Slabs')
pl.ylabel('Temperature Profile')

plt.legend(loc='best')
plt.savefig('Experiments.png')

plt.show()
