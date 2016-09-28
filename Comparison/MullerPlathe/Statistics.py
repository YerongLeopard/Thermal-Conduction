import numpy as np
import scipy.optimize as optimization
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
import scipy.optimize as optimization
def func(x,a,b):
  """linear function""" 
  return (a*x+b)
##Slab 0-10
KappaFile = open('KappaFile.dat', 'w')
ErrorFile = open('ErrorFile.dat','w')
temperature=pl.loadtxt('slab.dat')
error=pl.loadtxt('error.dat')
flux=pl.loadtxt('flux.dat')
i=0
x=np.arange(11)*30.1762008666992/20
x0=[0.0,0.0]##the initial guess of the parameters
#kappa=[]
#kappa_error=[]
for e in error[:][:]:
  y=temperature[:][i]
  f=flux[i]
  para,corv=optimization.curve_fit(func,x,y,x0,e)
  k=f/para[0]/(10.0586996078491**2)/2
  kerror=np.sqrt(corv[0][0])*k/para[0]
  print>> KappaFile,k
  print>> ErrorFile,kerror
#  kappa=kappa+[k]
#  kappa_error=kappa_error+[np.sqrt(corv[0][0])*k/para[0]]
  i=i+1
KappaFile.close()
ErrorFile.close()

##Slab 10-19
KappaFile = open('KappaFile2.dat','w')
ErrorFile = open('ErrorFile2.dat','w')
temperature=pl.loadtxt('slab2.dat')
error=pl.loadtxt('error2.dat')
i=0
x=np.arange(9)*30.1762008666992/20
x0=[0.0,0.0]##the initial guess of the parameters
#kappa=[]
#kappa_error=[]
for e in error[:][:]:
  y=temperature[:][i]
  f=flux[i]
  para,corv=optimization.curve_fit(func,x,y,x0,e)
  k=f/para[0]/(10.0586996078491**2)/2
  kerror=np.sqrt(corv[0][0])*k/para[0]
  #k=-k
  print>> KappaFile,-k
  print>> ErrorFile,kerror
  #kappa=kappa+[k]
  #kappa_error=kappa_error+[np.sqrt(corv[0][0])*k/para[0]]
  i=i+1
#print>> KappaFile,kappa
#print>> ErrorFile,kappa_error
KappaFile.close()
ErrorFile.close()

