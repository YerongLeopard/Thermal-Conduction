import numpy as np
import matplotlib;matplotlib.use("TkAgg");
import matplotlib.pyplot as plt
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
x=np.arange(27428-1200)
x=(x+1200)
for i in x:
  rerror1=rerror1+[error1[i]/data1[i]]
  rerror2=rerror2+[error2[i]/data2[i]]
  ratio=ratio+[(error1[i]*data2[i])/(error2[i]*data1[i])]
plt.figure('RelativeError')
# plt.plot(x*7,rerror1,'.r',label='origin',markersize=1)#old algorithm
plt.plot(x*7,rerror1,'.r',markersize=3)
# plt.plot(x*7,rerror2,'.b',label='new',markersize=1)#new algorithm
plt.plot(x*7,rerror2,'.b',markersize=3)
# plt.legend(loc='best') # No legend require
pl.xlabel('Step')
pl.ylabel('RelativeError')
plt.savefig('RelativeError')
plt.show()
