import scipy as sp;
import numpy as np;
import matplotlib;
try:
   from PyQt4 import QtCore, QtGui
except:
   matplotlib.use("TkAgg")
import pylab as pl;
def discrete_func(List, kappa, MinPath, scale = 1):     # discrete_func1 takes wavevctor q as the major input
  Debyek=0.985 # Debye cut-off wavelength depend merely on the density
  p=2; # p=2 has the real physics 
  N1=6
  N2=6
  ans=[]
  for N3 in List: 
    sum=0.0  
    i=0
    q1=2*sp.pi/N3  # minimum wavelength in this stucture 
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
          if ((-2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3>-3*N1*N2*N3)and(-2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3<3*N1*N2*N3+1)\
           and(+2*N2*N3*n1-2*N1*N3*n2+2*N1*N2*n3>-3*N1*N2*N3)and(+2*N2*N3*n1-2*N1*N3*n2+2*N1*N2*n3<3*N1*N2*N3+1)\
           and(+2*N2*N3*n1+2*N1*N3*n2-2*N1*N2*n3>-3*N1*N2*N3)and(+2*N2*N3*n1+2*N1*N3*n2-2*N1*N2*n3<3*N1*N2*N3+1)\
           and(+2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3>-3*N1*N2*N3)and(+2*N2*N3*n1+2*N1*N3*n2+2*N1*N2*n3<3*N1*N2*N3+1)):
             qc =q1*scale;  # calculating \kappa(qc=q1/4)
             x1=n1/float(N1)
             x2=n2/float(N2)
             x3=n3/float(N3)             
             sq=x1**2+x2**2+x3**2
             if (0!=sq):
               i+=1
               v=(n1/float(N1)+n2/float(N2)-n3/float(N3))/np.sqrt(3.0)
               k=np.sqrt(x1**2+x2**2+x3**2)
#               if (v>k):  # DEBUG
#                  print v,k # DEBUG
               result=result+((v**2)/sq)*((Debyek/k)**p)/(1+4*((np.sin(qc/2))**2)*(MinPath*(v/k)*((Debyek/k)**p)+(MinPath*(v/k)*((Debyek/k)**p))**2))
    result=result/(N3*N1*N2*4)*(3-p)*((np.cos(qc/2))**2)
    ans.append(result*kappa)
    # print i,'==',92*N2*4*N3,'(N=',N3,')' # DEBUG
  return ans
