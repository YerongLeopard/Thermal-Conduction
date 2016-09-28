import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
import scipy.optimize as optimization
def func(x,a,b): 
  return (a*np.exp(x*b))
data1= pl.loadtxt('KappaFile.dat')
error1=pl.loadtxt('ErrorFile.dat')
data2= pl.loadtxt('New.dat')
error2=pl.loadtxt('NewError.dat')
rerror1=[]
rerror2=[]
ratio=[]
x=np.arange(27428-1000)
x=(x+1000)
#x0=[0.0,0.0]
for i in x:
  rerror1=rerror1+[error1[i]/data1[i]]
  rerror2=rerror2+[error2[i]/data2[i]]
  ratio=ratio+[(error1[i]*data2[i])/(error2[i]*data1[i])]
plt.figure('RelativeError')
plt.plot(x*7,rerror1,'.r',label='origin',markersize=1)#old algorithm
plt.plot(x*7,rerror2,'.b',label='new',markersize=1)#new algorithm
plt.legend(loc='best')
pl.xlabel('Step')
pl.ylabel('RelativeError')
plt.savefig('RelativeError')

#plt.figure('ratio of relative error vs time')
#plt.plot(x,ratio,'.')
#plt.savefig('ratio of relative error vs time')

plt.figure('kappa vs time')
x=np.arange(27428-755)
plt.plot(x,data1[755:27428],'-')
plt.plot(x,data2[755:27428],'.')
plt.savefig('kappa vs time')
plt.show()
