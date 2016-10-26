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
   loop = np.arange(1,length,2) # Only do the summing up over odd nq's
   for nq in loop:
     q = q1 * nq
     if isInt:
        factor = sp.cos(q/2)
     divider = float(length)/float(nq)
     kappa = discrete_func([length],1,0.615140781881,divider)[0] # ratio
     rev_keff += factor * sp.sin(q*length/4)/sp.sin(q/2)/kappa
   rev_keff = 2*rev_keff/length
   keff = 1/rev_keff
   return keff
def main():
  L=list(np.arange(1500, 2000, 2 * 100));  # length has to be an even number 
  for l0 in L:
      File = open('k_eff_tmp.dat', 'a');
      keff=K_eff(l0)
      print >> File, l0, keff
      print l0,' ', keff
      File.close() 
  print("DEBUG")
  return


if __name__ == "__main__":
   DATA = np.loadtxt('6by6.dat') # Load data from the calculated results
main()
