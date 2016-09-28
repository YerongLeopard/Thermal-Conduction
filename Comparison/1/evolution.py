import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
import scipy.optimize as optimization
def func(x,a,b): 
  return (a*np.exp(x*b))
pl.rc('axes', linewidth=1.2)
fontsize=14
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')

ax.xaxis.set_tick_params(size=15,width=1.2)
ax.yaxis.set_tick_params(size=15,width=1.2)
data1= pl.loadtxt('KappaFile.dat')
error1=pl.loadtxt('ErrorFile.dat')
data2= pl.loadtxt('New.dat')
error2=pl.loadtxt('NewError.dat')
rerror1=[]
rerror2=[]
ratio=[]
x=np.arange(27428-1200)
x=(x+1200)
#x0=[0.0,0.0]
for i in x:
  rerror1=rerror1+[error1[i]/data1[i]]
  rerror2=rerror2+[error2[i]/data2[i]]
  ratio=ratio+[(error1[i]*data2[i])/(error2[i]*data1[i])]
#plt.figure('RelativeError')
plt.plot(x*(7.0)/50000,rerror1,'.r',markersize=3)#old algorithm
plt.plot(x*(7.0)/50000,rerror2,'.b',markersize=3)#new algorithm
#plt.legend(loc='best')
pl.xlabel('Step $\\times$ 50,000',fontsize=23)
pl.ylabel('Relative Standard Deviation',fontsize=23)
plt.xticks(range(0,5,1))
plt.yticks(np.arange(0,0.11,0.01))
plt.savefig('RelativeError')

#plt.figure('ratio of relative error vs time')
#plt.plot(x,ratio,'.')
#plt.savefig('ratio of relative error vs time')

#plt.figure('kappa vs time')
#x=np.arange(27428-755)
#plt.plot(x,data1[755:27428],'-')
#plt.plot(x,data2[755:27428],'.')
#plt.figure('kappa vs time')
plt.show()
