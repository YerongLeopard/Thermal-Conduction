import scipy as sp
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import mpl_toolkits.mplot3d.axes3d as p3
p=2
Debyek=0.985
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
      #       k1=n2/12.0+n3/float(N)/2#k vector in reciprocal space
      #       k2=n1/12.0+n3/float(N)/2
      #       k3=n1/12.0+n2/12.0
 #            bk1=n2*N+n3*6
 #            bk2=n1*N+n3*6
 #            bk3=n1*N+n2*N
             x1=n1/float(N1)
             x2=n2/float(N2)
             x3=n3/float(N3)             
 #            bx1=(bk2+bk3-bk1)
 #            bx2=(bk1+bk3-bk2)
 #            bx3=(bk1+bk2-bk3)
 #            check1=(-bx1+bx2+bx3)/3.0
 #            check2=(+bx1-bx2+bx3)/3.0
 #            check3=(+bx1+bx2-bx3)/3.0
        
             ##modify
             '''
             if check1>6*N:
                bk1-=12*N
                k1-=1
             if (check1<-6*N or check1==-6*N):
                bk1+=12*N
                k1+=1

             if check2>6*N:
                bk2-=12*N
                k2-=1
             if (check2<-6*N or check2==-6*N):
                bk2+=12*N
                k2+=1

             if check3>6*N:
                bk3-=12*N
                k3-=1
             if (check3<-6*N or check3==-6*N):
                bk3+=12*N
                k3+=1
             ##modify
             x1=(k2+k3-k1)
             x2=(k1+k3-k2)
             x3=(k1+k2-k3)

             bx1=(bk2+bk3-bk1)
             bx2=(bk1+bk3-bk2)
             bx3=(bk1+bk2-bk3)

             check1=(-bx1+bx2+bx3)/3.0
             check2=(+bx1-bx2+bx3)/3.0
             check3=(+bx1+bx2-bx3)/3.0

#             if i==861:
#               print '\n',i 
#               print 'k',k1,k2,k3
#               print 'x',x1,x2,x3
#               print 'check',check1/12.0/N,check2/12.0/N,check3/12.0/N
             
             if check1>6*N:
                print 'error'

             if (check1<-6*N or check1==-6*N):
                print 'error'

             if check2>6*N:
                print 'error'

             if (check2<-6*N or check2==-6*N):
                print 'error'

             if check3>6*N:
                print 'error'

             if (check3<-6*N or check3==-6*N):
                print 'error'
             
#             bx1=(bk2+bk3-bk1)
#             bx2=(bk1+bk3-bk2)
#             bx3=(bk1+bk2-bk3) 
             '''
#             ax.scatter(x1,x2,x3,s=10)
             sq=x1**2+x2**2+x3**2
#             bsq=bx1**2+bx2**2+bx3**2
             if (0!=sq):
               v=(n1/float(N1)+n2/float(N2)-n3/float(N3))/np.sqrt(3.0)
               k=np.sqrt(x1**2+x2**2+x3**2)
#               bk=np.sqrt(bx1**2+bx2**2+bx3**2)
 #              if (k<v):
 #                 print k
     #          print k*12*N,bk
               result=result+((v**2)/sq)*((Debyek/k)**p)/(1+4*((np.sin(q0/2))**2)*(MinPath*(v/k)*((Debyek/k)**p)+(MinPath*(v/k)*((Debyek/k)**p))**2))
            #   print ((v**2)/sq)*((Debyek/k)**p)*((np.cos(q0/2))**2)/(1+4*((np.sin(q0/2))**2)*(MinPath*(v/k)*((Debyek/k)**p)+(MinPath*(v/k)*((Debyek/k)**p))**2)),'='
               #print '$v/k$==',v/k,sum
    result=result/(N3*N1*N2*4)*(3-p)*((np.cos(q0/2))**2)
    ans.append(result*kappa)
#  ax.set_xlabel('X')
#  ax.set_ylabel('Y')
#  ax.set_zlabel('Z')
#  fig.add_axes(ax)
    print(i,'==',N1*N2*4*N3,'(N=',N3,')')
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

w0=0.5
w=w0
data=pl.loadtxt('data-all.dat')
q=[]
kappa=[]
error=[]
x0=np.array([20,1.75])
for l in data:
  q.append((2*sp.pi/l[0])**w0)
  kappa.append(l[1])
  error.append(l[2])
q=np.array(q)
kappa=np.array(kappa)
error=np.array(error)
#plt.plot(q,kappa,'b.',markersize=16)
#plt.errorbar(q,kappa,yerr=(error),linestyle="None",linewidth=3)

line=pl.loadtxt('6by6.dat')
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
plt.plot(x,y,'r-',linewidth=2,label='$6\\times 6$')


line=pl.loadtxt('4by4.dat')
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
plt.plot(x,y,'b-',linewidth=2,label='$4 \\times 4$')

line=pl.loadtxt('2by2.dat')
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
plt.plot(x,y,'g-',linewidth=2,label='$2 \\times 2$')

pl.legend(loc='best',prop={'size':20})
#pl.title("6s by 6")
pl.xlabel('$\sqrt{qa}$',fontsize=25)
pl.ylabel('$\kappa(q)/\kappa(0)$',fontsize=30)
#pl.xlabel('$\sqrt{qa}%f$'%w0,fontsize=25)
#pl.ylabel('$\kappa(q)$',fontsize=30)
pl.xticks(fontsize=14)
pl.yticks(fontsize=14)
print('finished \n')
plt.show()
