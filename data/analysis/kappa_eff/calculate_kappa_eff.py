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
   # print isInt # DEBUG
   # Nq = np.arange(1, length, 2) # Only do the summing up over odd nqs
   Nq = [item for item in np.arange(-int(length/2)+1 , int(length/2), 1) if item % 2]
   # print(Nq, type(Nq)) # DEBUG 
   kappa = discrete_func([length],1,0.615140781881,Nq)[0]# ratio
   # print kappa, 'kappa' # DEBUG
   if isInt:
     tmp = [sp.sin(nq*sp.pi/2)*sp.cos(nq*sp.pi/length)/sp.sin(nq*sp.pi/length)/k for nq,k in zip(Nq, kappa)]
   else:
     tmp = [sp.sin(nq*sp.pi/2)*1                      /sp.sin(nq*sp.pi/length)/k for nq,k in zip(Nq, kappa)]
   # print tmp[0: 5], 'tmp' # DEBUG
   # print tmp[int(length/4):], 'tmp'  # DEBUG
   # print [1/sp.sin(nq*sp.pi/length) for nq in Nq][0: int(length/4)], '1/sin' # DEBUG
   # print [sp.sin(nq*sp.pi/2) for nq in Nq][0: int(length/4)], 'sp.sin(nq*sp.pi/2)' # DEBUG
   # print [sp.sin(nq*sp.pi/2) for nq in Nq][int(length/4):], 'sp.sin(nq*sp.pi/2)' # DEBUG
   # list1 = [sp.sin(nq*sp.pi/2) for nq in Nq] # DEBUG
   # list2 = [1/sp.sin(nq*sp.pi/length) for nq in Nq] # DEBUG
   # list3 = [sp.cos(nq*sp.pi/length) for nq in Nq] # DEBUG
   # print list3[0:5] # DEBUG
   # print "test1", list1[0]*list2[0]*list3[0]/kappa[0], list1[1] * list2[1]*list3[1]/kappa[1] # DEBUG
   rev_keff = 2/float(length)*np.sum(tmp)
   keff = 1/rev_keff
   return keff
def main():
  L=list(np.arange(2600, 3000 +1, 200));  # length has to be an even number 
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
