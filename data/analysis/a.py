import pylab as pl
import scipy.optimize as optimization
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy as sp
import math
from scipy.integrate import dblquad,quad
Debyek=0.985
p=2
N1=6
N2=6
def continuous_func(q,kappa,MinPath):
    ans=[] 
    for q0 in q:
       a,b= dblquad(lambda x,theta:((np.cos(theta))**2)*np.sin(theta)*(x**(-p+2))/(1+4*(np.sin(q0/2))**2*(MinPath*(x**(-p))*np.cos(theta)+(MinPath*(x**(-p))*np.cos(theta))**2)),0,sp.pi,lambda x:0,lambda x:1)
       ans.append(a*((np.cos(q0/2))**2)*3*((3-p)/2.0)*kappa)
    return  ans
#  return answer
#--------------------------------------------------------------------------------------------------------------------------

def discrete_func1(q,kappa,MinPath):
#  fig=pl.figure()
#  ax = p3.Axes3D(fig)
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
    #k=int(N/2)
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
             ax.scatter(x1,x2,x3,s=10)
             sq=x1**2+x2**2+x3**2
             if (0!=sq):
               v=(n1/float(N1)+n2/float(N2)-n3/float(N3))/np.sqrt(3.0)
               k=np.sqrt(x1**2+x2**2+x3**2)
               result=result+((v**2)/sq)*((Debyek/k)**p)/(1+4*((np.sin(q0/2))**2)*(MinPath*(v/k)*((Debyek/k)**p)+(MinPath*(v/k)*((Debyek/k)**p))**2))
    result=result/(N3*N1*N2*4)*(3-p)*((np.cos(q0/2))**2)
    ans.append(result*kappa)
#  ax.set_xlabel('X')
#  ax.set_ylabel('Y')
#  ax.set_zlabel('Z')
#  fig.add_axes(ax)
    print i,'==',N1*N2*4*N3,'(N=',N3,')'
  return ans
#--------------------------------------------------------------------------------------------------------------------------


def discrete_func2(q,L1,kappa,MinPath):
#  fig=pl.figure()
#  ax = p3.Axes3D(fig)

  ans=[]
  for N1 in L1:
    N2=N1
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
    print i,'==',N1*N2*4*N3,'(N3=',N3,')','N1',N1
  return ans
pl.rc('axes', linewidth=2)
fontsize = 14
ax = pl.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
ax.xaxis.set_tick_params(size=7,width=2)
ax.yaxis.set_tick_params(size=7,width=2)
#print discrete_func2([2*np.pi/10000],[60],1,-0.525031555126)
'''
############################   continuous:Debye 
q=np.arange(0,0.9,0.01)
p=0
plt.plot(q,continuous_func(q*np.pi,1,4,10),'-b',label="p="+str(p))
p=1
plt.plot(q,continuous_func(q*np.pi,1,4,10),'--g',label="p="+str(p))
p=2
plt.plot(q,continuous_func(q*np.pi,1,4,10),'.r',label="p="+str(p))
plt.legend(loc='best')
pl.xlabel('$q$',fontsize=20)
pl.ylabel('$\kappa(q)$/$\kappa_D$',fontsize=20)
pl.show()
'''
'''
##########################    continuous log:Debye
pl.figure(1)
q=np.arange(0,0.9,0.001)
p=0
plt.plot(q,np.log10(continuous_func(q*np.pi,1,10)),'.b',label="p="+str(p))
p=1
plt.plot(q,np.log10(continuous_func(q*np.pi,1,10)),'--g',label="p="+str(p))
p=2
plt.plot(q,np.log10(continuous_func(q*np.pi,1,10)),'-r',label="p="+str(p))
pl.ylim([-3,0])
plt.legend(loc='best')
pl.xlabel('$qd/\pi$',fontsize=20)
pl.ylabel('$log(\kappa(q)$/$\kappa_D$)',fontsize=20)
pl.show()



##########################     discrete :BZ
pl.figure(2)
q=np.arange(0.008,0.9,0.001)
p=0
plt.plot(q,np.log10(discrete_func1(q*np.pi,1,10)),'.b',label="p="+str(p))
p=1
plt.plot(q,np.log10(discrete_func1(q*np.pi,1,10)),'--g',label="p="+str(p))
p=2
plt.plot(q,np.log10(discrete_func1(q*np.pi,1,10)),'-r',label="p="+str(p))

pl.ylim([-3,0])
plt.legend(loc='best')
pl.xlabel('$qd/\pi$',fontsize=20)
`	pl.ylabel('$log(\kappa(q)$/$\kappa_D)$',fontsize=20)
pl.show()
'''

'''
##########################     discrete :BZ Small D
q=np.arange(0,0.9,0.01)
p=0
plt.plot(q,np.log10(discrete_func1(q*np.pi,1,0.4,1.75)),'-b',label="p="+str(p))
p=1
plt.plot(q,np.log10(discrete_func1(q*np.pi,1,0.4,1.75)),'--g',label="p="+str(p))
p=2
plt.plot(q,np.log10(discrete_func1(q*np.pi,1,0.4,1.75)),'.r',label="p="+str(p))
pl.ylim([-2,0])
plt.legend(loc='best')
pl.xlabel('$qd/\pi$',fontsize=20)
pl.ylabel('$log(\kappa(q)$/$\kappa_D)$',fontsize=20)
pl.show()
'''
#####################plot cross section
#x = np.linspace(0, 2 * np.pi, 400)
#y = np.sin(x ** 2)
#f, axarr = plt.subplots(2, sharex=True)
#axarr[0].plot(x, y)
#axarr[0].set_title('Sharing X axis')
#axarr[1].scatter(x, y)
'''
ax = pl.subplot(1,1,1)
p1, = ax.plot([1,2,3], label="line 1")
p2, = ax.plot([3,2,1], label="line 2")
p3, = ax.plot([2,3,1], label="line 3")

handles, labels = ax.get_legend_handles_labels()

# reverse the order
ax.legend(handles[::-1], labels[::-1])

# or sort them by labels
#import operator
#hl = sorted(zip(handles, labels),
            key=operator.itemgetter(1))
#handles2, labels2 = zip(*hl)

ax.legend(handles2, labels2)
'''
pl.figure(1)
#f, axarr = plt.subplots(2, sharex=True)
############################Writing the data in:cross section
'''
File = open('cross.dat', 'w')
n=2000
number=7
refer=discrete_func2([2*np.pi/n],[6*(2**(number-1))],1,10)
l=[]
x=[]
y=[]
for w in range(0,(number-1),1):
  element=6*(2**w)
  l.append(element)
  x.append(np.log10(element))
y=discrete_func2([2*np.pi/n],l,1,10)
for w in range(len(y)):
  print>> File,n,6*(2**w),y[w]  
print>> File,n,6*(2**(number-1)),refer[0]
File.close()
'''

########################plotting:cross section
x=[]
y=[]
data=pl.loadtxt('cross1.dat')
for element in data[0:]:
   y.append(element[2]-0.01)
   x.append(math.log(element[1],5))
refer=data[-1][2]-0.01
#y.append(refer)
#x.append(np.log10(data[-1][1]))
for w in range(len(y)):
  y[w]=y[w]/refer
math
plt.plot(x,y,'.b',markersize=15)
pl.xlabel('Width of Cross Section $N_1$',fontsize=20)
pl.ylabel('$\kappa_{N_1,N_3}/\kappa_{max}$',fontsize=30)
plt.xticks([i for i in range(1,4,1)],['%i'%(5**i) for i in range(1,4,1)])

plt.show()


'''
##########################    discrete log:Debye

q=np.arange(0,0.9,0.005)

p=0
plt.plot(q,np.log10(discrete_func2(q*np.pi,1,2,10)),'-b',label="p="+str(p))
p=1
plt.plot(q,np.log10(discrete_func2(q*np.pi,1,2,10)),'--g',label="p="+str(p))
p=2
plt.plot(q,np.log10(discrete_func2(q*np.pi,1,2,10)),'.r',label="p="+str(p))

plt.legend(loc='best')
pl.xlabel('$qd/\pi$',fontsize=20)
pl.ylim([-3,0])
pl.ylabel('$log(\kappa(q)$/$\kappa_D$)',fontsize=20)
pl.show()'''
