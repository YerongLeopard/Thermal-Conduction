import scipy as sp
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import mpl_toolkits.mplot3d.axes3d as p3
p=2
Debyek=0.985
def discrete_func1(q,kappa,MinPath):
  N1=6
  N2=6
  ans=[]
  for q0 in q: 
    sum=0.0  
    i=0
    N3=6
    if (abs(q0)>10e-10):
      N3=round(2*sp.pi/q0,0)
      N3=int(N3)
    first1=-N1+1
    last1=N1+1
    first2=-N2+1
    last2=N2+1
    first3=-N3+1
    last3=N3+1
    result=0.0
    for n1 in range(first1,last1):
      for n2 in range(first2,last2):
        for n3 in range(first3,last3):
          if ((-2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3>-3*N1*N2*N3)and(-2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3<3*N1*N2*N3+1)and(+2*N2*N3*n1-2*N1*N3*n2+2*N1*N2*n3>-3*N1*N2*N3)and(+2*N2*N3*n1-2*N1*N3*n2+2*N1*N2*n3<3*N1*N2*N3+1)and(+2*N2*N3*n1+2*N1*N3*n2-2*N1*N2*n3>-3*N1*N2*N3)and(+2*N2*N3*n1+2*N1*N3*n2-2*N1*N2*n3<3*N1*N2*N3+1)and(+2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3>-3*N1*N2*N3)and(+2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3<3*N1*N2*N3+1)):
             i+=1
             x1=n1/float(N1)
             x2=n2/float(N2)
             x3=n3/float(N3)             
             sq=x1**2+x2**2+x3**2
             if (0!=sq):
               v=(n1/float(N1)+n2/float(N2)-n3/float(N3))/np.sqrt(3.0)
               k=np.sqrt(x1**2+x2**2+x3**2)
               result=result+((v**2)/sq)*((Debyek/k)**p)/(1+4*((np.sin(q0/2))**2)*(MinPath*(v/k)*((Debyek/k)**p)+(MinPath*(v/k)*((Debyek/k)**p))**2))
    result=result/(N3*N1*N2*4)*(3-p)*((np.cos(q0/2))**2)
    ans.append(result*kappa)
    print (i,'==',N1*N2*4*N3-1,'(N=',N3,')')
  return ans

pl.rc('axes', linewidth=1.2)
fontsize =10
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
ax.xaxis.set_tick_params(size=15,width=1.2)
ax.yaxis.set_tick_params(size=15,width=1.2)
line=pl.loadtxt('6by6real.dat')
w0=0.5
w=w0
data=pl.loadtxt('data-all.dat')
q=[]
kappa=[]
error=[]
x0=np.array([23,0.61514])
for l in data:
  q.append((2*sp.pi/l[0]))
  kappa.append(l[1])
  error.append(l[2])
q=np.array(q)
kappa=np.array(kappa)
error=np.array(error)
para = optimization.curve_fit(discrete_func1, q, kappa, x0, error)
print para, " para"
#plt.plot(q,kappa,'b.',markersize=16)
#plt.errorbar(q,kappa,yerr=(error),linestyle="None",linewidth=3)

q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append(np.sqrt(2*sp.pi/l[0])*(l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0)
plt.plot(x,y,'r-',linewidth=2,label='$6\\times 6$')
pl.legend(loc='best',prop={'size':20})
pl.xlabel('$\sqrt{qa}$',fontsize=25)
pl.ylabel('$\kappa(q)/\kappa(0)\sqrt{qa}$',fontsize=30)
pl.xticks(fontsize=14)
pl.yticks(fontsize=14)
plt.show()
