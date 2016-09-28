import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
data= pl.loadtxt('Comparison.dat')
plt.figure(1)
x=[]
y=[]
y1=[]
y2=[]
for i in range(0,31):
  x=x+[data[i][0]]
  y=y+[data[i][1]]
  y1=y1+[data[i][2]]
  y2=y2+[data[i][3]]
#plt.title='The ratio of the time step of new algorithm over the old one'  
plt.scatter(x,y,s=20)
plt.plot([0],[1],'.r',markersize=10)
pl.xticks([0.000,0.010,0.020,0.030,0.040,0.050,0.060])
pl.ylabel('Ratio of time steps')
pl.xlabel('Relative error')
plt.savefig('Comparison.png')
plt.show()
