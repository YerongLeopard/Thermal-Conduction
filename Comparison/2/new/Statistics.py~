import numpy as np
import scipy.optimize as optimization
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
import scipy.optimize as optimization
import math as m
PI=3.1415926535898
amplitude=1.265*4
scale=30.1762008666992/20
KappaFile = open('New.dat', 'w')
ErrorFile = open('NewError.dat','w')
def func(x,a,b): 
 return (-a*np.cos(2*PI*(x+1)/20)+b) 
#plotting
temperature= pl.loadtxt('slab.dat')
x =np.arange(20)
error=pl.loadtxt('error.dat')
x0=np.array([0.0,0.0])
i=0
for e in error[:][:]:
  y=temperature[:][i]
  para,corv=optimization.curve_fit(func,x,y,x0,e)
  T=para[0]
  k=amplitude*scale/T/4/(np.sin(PI/20)**2)/(10.0586996078491**2)
  Terror=np.sqrt(corv[0][0])
  kerror=Terror/T*k
  print>> KappaFile,k
  print>> ErrorFile,kerror
  i=i+1
KappaFile.close()
ErrorFile.close()
