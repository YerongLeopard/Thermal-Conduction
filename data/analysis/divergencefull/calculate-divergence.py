from analysis import *
import scipy as sp;
import numpy as np;
import matplotlib;
try:
   from PyQt4 import QtCore, QtGui
except:
   matplotlib.use("TkAgg")
import pylab as pl;
import matplotlib.pyplot as plt;
import scipy.optimize as optimization;
import mpl_toolkits.mplot3d.axes3d as p3;
L=list(np.arange(20,21, 1));
for l0 in L:
  File = open('6by.dat', 'a');
  y=discrete_func([l0],1,0.615140781881,[1])[0] # ratio
  print >> File,l0,y[0]
  print l0,' ', y[0]
  File.close()
print '\nfinished.\n'

