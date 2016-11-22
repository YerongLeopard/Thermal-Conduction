import numpy as np
import matplotlib;
try:
   from PyQt4 import QtCore, QtGui
except:
   matplotlib.use("TkAgg")
import pylab as pl
import matplotlib.pyplot as plt
import scipy.optimize as optimization
# Customize frame
pl.rc('axes', linewidth=2)
data1= pl.loadtxt('KappaFile.dat')
error1=pl.loadtxt('ErrorFile.dat')
data2= pl.loadtxt('New.dat')
error2=pl.loadtxt('NewError.dat')
rerror1=[]
rerror2=[]
ratio=[]
x=np.arange(27428-1200) # Original 
# x = np.arange(1000) # DEBUG
x=(x+1200)
x_position = [25000, 50000, 75000, 100000, 125000, 150000, 175000, 200000]
x_label = ['', '1','', '2', '', '3', '', '4']
for i in x:
  rerror1=rerror1+[error1[i]/data1[i]]
  rerror2=rerror2+[error2[i]/data2[i]]
  ratio=ratio+[(error1[i]*data2[i])/(error2[i]*data1[i])]
plt.figure('RelativeError')
# plt.plot(x*7,rerror1,'.r',label='origin',markersize=1)#old algorithm
plt.plot(x*7,rerror1,'.r',markersize=5)
# plt.plot(x*7,rerror2,'.b',label='new',markersize=1)#new algorithm
plt.plot(x*7,rerror2,'.b',markersize=5)
# plt.legend(loc='best') # No legend require

ax = plt.gca()
# Customize tick marks
major_ticks = np.arange(0, 200001, 50000)
minor_ticks = np.arange(0, 200001, 25000)
ax.tick_params(width=2, length=8) # setting linth and width of tick marks

# Customize tick labels
for idx, tick in enumerate(ax.xaxis.get_major_ticks()):
    tick.label1.set_fontsize(20)
    tick.label1.set_fontweight('bold')

for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
    tick.label1.set_fontweight('bold')



pl.xlim(0 , 200000)
plt.xticks(x_position, x_label)
pl.ylim(0.00 , 0.10)
pl.xlabel('Step($ \\times 50,000$)', fontsize=25)
pl.ylabel('Relative Standard Deviation', fontsize=25)
plt.savefig('RelativeError')
plt.show()
