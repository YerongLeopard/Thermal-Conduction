import sys;
sys.path.append("../include/")
from analysis import *
import numpy as np
import scipy as sp
import matplotlib
try:
   from PyQt4 import QtCore, QtGui
except:
   matplotlib.use("TkAgg")
def K_eff(length):
   isInteger = False; factor = 1.0
   q1 = 2 * sp.pi / length 
   rev_keff = 0
   if (0 != length%2):
     print("length has to be an even number") # Error message 
     quit()
   elif (0 != length%4):  # Ns/4 is a half integer
     isInt = False
   else:                  # Ns/4 is an integer
     isInt = True
   Nq = [item for item in np.arange(-int(length/2)+1 , int(length/2), 1) if item % 2]
   kappa = discrete_func([length],1,0.615140781881,Nq)[0]# ratio
   if isInt:
     tmp = [sp.sin(nq*sp.pi/2)*sp.cos(nq*sp.pi/length)/sp.sin(nq*sp.pi/length)/k for nq,k in zip(Nq, kappa)]
   else:
     tmp = [sp.sin(nq*sp.pi/2)*1                      /sp.sin(nq*sp.pi/length)/k for nq,k in zip(Nq, kappa)]
   rev_keff = 2/float(length)*np.sum(tmp)
   keff = 1/rev_keff
   return keff
def main():
  L=list(np.arange(100, 2800 +1, 100));  # length has to be an even number 
  for l0 in L:
      File = open('k_eff_tmp4.dat', 'a');
      keff=K_eff(l0)
      print >> File, l0, keff
      print l0,' ', keff
      File.close() 
  print("DEBUG")
  return


if __name__ == "__main__":
   DATA = np.loadtxt('6by6.dat') # Load data from the calculated results
main()
