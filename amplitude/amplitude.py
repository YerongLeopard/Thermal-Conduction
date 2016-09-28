import numpy as np
import scipy.optimize as optimization
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
def func(x,a,b): 
  return (a*x+b)
pl.rc('axes', linewidth=1.2)
fontsize = 12
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
ax.xaxis.set_tick_params(size=7,width=1.2)
ax.yaxis.set_tick_params(size=7,width=1.2)
#plotting1
plt.figure(1)
x =np.arange(11)
y= pl.loadtxt('slab.dat')
flux=pl.loadtxt('flux.dat')
err=pl.loadtxt('error.dat')
y= y[-1]
flux=flux[-1]
err=err[-1]
x0=np.array([0.0,0.0])
para,corv=optimization.curve_fit(func,x,y,x0,err)
k=para[0]
b=para[1]
#calculation
x1=x*(30.1762008666992/20)
para,corv=optimization.curve_fit(func,x1,y,x0,err)
k_real=para[0]
kappa=flux/k_real/(10.0586996078491**2)/2
kappa_devi=np.sqrt(corv[0][0])*kappa/k_real
#calculation
plt.plot(x,y,'.b',markersize=10)
plt.plot( x,k*x+b,'-r',linewidth=2,label=('$\kappa=$'+str(round(kappa,2))+'('+str(0.14)+')'))
plt.errorbar(x,y,yerr=(err),linestyle="None")


#plotting2
x =np.arange(10)
x=x+10
y= pl.loadtxt('slab2.dat')
flux=pl.loadtxt('flux.dat')
err=pl.loadtxt('error2.dat')
y= y[-1]
flux=flux[-1]
err=err[-1]
x0=np.array([0.0,0.0])
para,corv=optimization.curve_fit(func,x,y,x0,err)
k=para[0]
b=para[1]
#calculation
x1=x*(30.1762008666992/20)
para,corv=optimization.curve_fit(func,x1,y,x0,err)
k_real=para[0]
kappa=flux/k_real/(10.0586996078491**2)/2
kappa_devi=np.sqrt(corv[0][0])*kappa/k_real
#calculation
plt.plot(x,y,'.b',markersize=10)
x=list(x)
x.append(20)
x=np.array(x)
plt.plot(x,k*x+b,'-r',linewidth=2)#,label=('Original 2nd Fit:$\kappa=$'+str(round(-kappa,2))+'('+str(round(kappa_devi,2))+')'))

#plotting3
PI=3.1415926535898
e=1.265
def func(x,a,b): 
 return (-a*np.cos(2*PI*(x)/20)+b) 
x =np.arange(21)
y= pl.loadtxt('slabnew.dat')
err=pl.loadtxt('errornew.dat')

y=y[-1]
err=err[-1]
y=list(y)
err=list(err)

tmp=[y[-1]]
for i in range(len(y)):
  tmp.append(y[i])
y.append(1)
for i in range(len(y)):
  y[i]=tmp[i]

tmp=[err[-1]]
for i in range(len(err)):
  tmp.append(err[i])
err.append(1)
for i in range(len(err)):
  err[i]=tmp[i]

x0=np.array([0.0,0.0])
para,corv=optimization.curve_fit(func,x,y,x0,err)
T=para[0]
T0=para[1]
#calculation
scale=30.1762008666992/20
kappa=e*scale/T/4/(np.sin(PI/20)**2)/(10.0586996078491**2)
Terror=np.sqrt(corv[0][0])
kerror=Terror/T*kappa
plt.plot(x,y,'.g',markersize=10)
#calculation
#plotting
x1=np.arange(20000)
x1=x1/1000.0
plt.plot(x1,func(x1,T,T0),'-b',linewidth=2,label='$\kappa=$'+'7.04(0.07)')
plt.errorbar(x,y,yerr=(err),linestyle="None")
x=np.arange(20)
#plt.xticks([i in range(20)],[i for i in range(20)])
plt.xticks([i for i in range(21)],[np.mod(i-10,20) for i in range(21)])
#plt.xticks([i in range(4)],['%i'%(10**i) for i in range(4)])
pl.legend(loc='best')
pl.xlabel('Slabs',fontsize=23)
pl.ylabel('Temperature',fontsize=23)

plt.show()
