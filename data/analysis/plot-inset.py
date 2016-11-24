import pylab as pl
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import scipy.optimize as optimization
fig = plt.figure()
pl.rc('axes', linewidth=1.2)
text10='10by10.dat'
text9='9by9.dat'
text8='8by8.dat'
text7='7by7.dat'
text6='6by6real.dat'
text5='5by5.dat'
text4='4by4.dat'
text3='3by3.dat'
text2='2by2.dat'
start=14
end=272
xx= np.linspace(0, 0.8, 10, endpoint=True) # extrapolation range
def func(x,a,b,c): 
  return (a+b*x+c*x**2) 
 
def func1(x,a,b,c): 
  return (a/x+b+c*x) 

def lin(x,b,c):
  return b+c*x
##################################
def func2(x,a,b,c,d,e,f): 
 return (a+b*x+c*x**2+d*x**3+e*x**4+f*x**5) 

def func3(x,a,b,c,d,e,f): 
 return (a/x+b+c*x+d*x**2+e*x**3+f*x**4) 
axes1 = fig.add_axes([0.12, 0.14, 0.8, 0.8],xlim=[0.0,0.8]) # main axes
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 20)
axes1.xaxis.set_tick_params(size=14,width=1.2)
axes1.yaxis.set_tick_params(size=14,width=1.2)
axes2 = fig.add_axes([0.35, 0.5, 0.4, 0.35],xticks=np.arange(0.0,0.2+0.05,0.05),xlim=[0.0,0.2],ylim=[0.0,0.6]) # inset axes
plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)
# main figure
w=0.5
line=pl.loadtxt(text10)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes1.plot(x,y,'m-',linewidth=2,label='$10\\times 10$')

line=pl.loadtxt(text9)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes1.plot(x,y,'c-',linewidth=2,label='$9\\times 9$')

line=pl.loadtxt(text8)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes1.plot(x,y,'c-',linewidth=2,label='$8\\times 8$')

line=pl.loadtxt(text7)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes1.plot(x,y,'c-',linewidth=2,label='$7\\times 7$')

line=pl.loadtxt(text6)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y] # temporary correction
axes1.plot(x,y,'r-',linewidth=2,label='$6\\times 6$')
'''
line=pl.loadtxt(text5)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y]
axes1.plot(x,y,'c-',linewidth=2,label='$5\\times 5$')
'''
line=pl.loadtxt(text4)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y]
axes1.plot(x,y,'b-',linewidth=2,label='$4\\times 4$')

line=pl.loadtxt(text3)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y]
axes1.plot(x,y,'k-',linewidth=2,label='$3\\times 3$')

line=pl.loadtxt(text2)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y]
axes1.plot(x,y,'g-',linewidth=2,label='$2\\times 2$')
#axes1.legend(loc='center left')
axes1.set_xlabel('$\sqrt{qa}$',fontsize=30)
axes1.set_ylabel('$\kappa(q)/\kappa_{0}$',fontsize=25)
#axes1.set_title('title')

# insert ##########################################################################################
num=10
w=0.5
line=pl.loadtxt(text10)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes2.plot(x,y,'m-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
#axes1.plot(x,func1(x,para[0],para[1],para[2]),'k--',markersize=1)
num-=1

line=pl.loadtxt(text9)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes2.plot(x,y,'m-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
#axes1.plot(x,func1(x,para[0],para[1],para[2]),'k--',markersize=1)
num-=1

line=pl.loadtxt(text8)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes2.plot(x,y,'c-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
#axes1.plot(x,func1(x,para[0],para[1],para[2]),'k--',markersize=1)
num-=1

line=pl.loadtxt(text7)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
#axes2.plot(x,y,'m-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
#axes1.plot(x,func1(x,para[0],para[1],para[2]),'k--',markersize=1)
num-=1

line=pl.loadtxt(text6)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y] # temporary DEBUG
axes2.plot(x,y,'r-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
x1=np.arange(0.025,0.79,0.025)
axes1.plot(x1,func1(x1,para[0],para[1],para[2]),'ro',linewidth=2)
print xx
print func1(xx,0,para[1],para[2])
axes1.plot(xx,lin(xx,para[1],para[2]),'r--',linewidth=2)
num-=1

line=pl.loadtxt(text5)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y] # temporary DEBUG
#axes2.plot(x,y,'c-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
#axes1.plot(x,func1(x,para[0],para[1],para[2]),'k--',linewidth=2)
#axes1.plot(x,func1(x,para[0],para[1],para[2],para[3]),'k--',markersize=1)
num-=1

line=pl.loadtxt(text4)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y] # temporary DEBUG
axes2.plot(x,y,'b-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
x1=np.arange(0.025,0.79,0.025)
axes1.plot(x1,func1(x1,para[0],para[1],para[2]),'bo',linewidth=2)
axes1.plot(xx,lin(xx,para[1],para[2]),'b--',linewidth=2)
#axes1.plot(x,func1(x,para[0],para[1],para[2],para[3]),'k--',markersize=1)
num-=1

line=pl.loadtxt(text3)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y] # temporary DEBUG
axes2.plot(x,y,'k-',linewidth=2)

x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
x1=np.arange(0.025,0.79,0.025)
axes1.plot(x1,func1(x1,para[0],para[1],para[2]),'ko',linewidth=2)
axes1.plot(xx,lin(xx,para[1],para[2]),'k--',linewidth=2)
#axes1.plot(x,func1(x,para[0],para[1],para[2],para[3]),'k--',markersize=4)
num-=1

line=pl.loadtxt(text2)
q=[]
y=[]
for l in line[:]:
  q.append(2*sp.pi/l[0])
  y.append((l[1])*np.sqrt(2*sp.pi/l[0]))#30.4962234594))
x=[]
for q0 in q:
     x.append(q0**w)
y=[element * 3 for element in y] # temporary DEBUG
axes2.plot(x,y,'g-',linewidth=2)
x=np.array(x)
y=np.array(y)
x0=np.array([0.0,0.0,0.0])
#x0=np.array([0.0,0.0,0.0,0.0])
x1=x[start:end]
y1=y[start:end]
para,corv=optimization.curve_fit(func,x1,y1,x0)
print num,para[0],para[1],para[2]
x1=np.arange(0.025,0.79,0.025)
axes1.plot(x1,func1(x1,para[0],para[1],para[2]),'go',linewidth=2)
axes1.plot(xx,lin(xx,para[1],para[2]),'g--',linewidth=2)

axes2.grid(True)
axes2.set_xlabel('$\sqrt{qa}$',fontsize=20)
axes2.set_ylabel('$\kappa(q)\sqrt{qa}/\kappa_{0}$',fontsize=20)
#axes1.plot(x,func1(x,para[0],para[1],para[2],para[3]),'k.',markersize=4)
#axes1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
axes1.legend(loc='best',fontsize=15)
plt.show()
