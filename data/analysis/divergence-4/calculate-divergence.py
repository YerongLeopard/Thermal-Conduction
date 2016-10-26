import sys;
sys.path.append("../include/")
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
L=list(np.arange(400,12000*4, 4));
for l0 in L:
  File = open('6by.dat', 'a');
  y=discrete_func([l0],1,0.615140781881,2) # ratio
  print >> File,l0,y[0]
  print l0,' ', y[0]
  File.close()
print '\nfinished. \nkappa(q/4)\n'
