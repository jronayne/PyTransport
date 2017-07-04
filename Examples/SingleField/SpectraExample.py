####################################### MTeasyPyStep simple example of basic functions ###########################################
from matplotlib import pyplot as plt
import time
import imp  
from pylab import *
import numpy as np
from scipy import interpolate 
import sympy as sym
import subprocess
############################################################################################################################################




############################################################################################################################################

#This file contains simple examples of using the PyTrans package for the single field example of Chen et al.
#It assumes the StepExampleSeptup file has been run to install a Step version of MTeasyPyStep
#It is recommended you restart the kernel to insure any updates to MTeasyPyStep are imported 

location = "/Users/David/Dropbox/MTeasyDist/MTeasy/" # this should be the location of the MTeasy folder 
sys.path.append(location) # sets up python path to give access to MTeasySetup

import MTeasySetup
MTeasySetup.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyStep as MTSE;  # import module  
import MTeasyScripts as MTS;
###########################################################################################################################################



# Example 
tols = np.array([10**-8,10**-8])

#################################### set initial field values and parameters  #############################################################

nF=MTSE.nF() # gets number of fields (useful check)
nP=MTSE.nP() # gets number of parameters needed (useful check)

fields = np.array([17.0])

params = np.zeros(nP)
params[0]=pow(10.,-5); params[1]=0.0018; params[2]=14.84; params[3]=0.022;

V = MTSE.V(fields,params) # calculate potential from some initial conditions
dV=MTSE.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)

initial = np.array([fields,-dV/np.sqrt(3.*V)]) # set initial conditions to be in slow roll
############################################################################################################################################




################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 40.0
t=np.linspace(Nstart, Nend, 1000)
back = MTSE.backEvolve(t, initial, params,tols,False)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g') # always good to inspect this plot to make sure its senisble
###########################################################################################################################################




################################## Calculate power spectrum and bispectrum for equilateral configs ########################################
fnlOut=np.array([])
kOut=np.array([])
PhiOut=np.array([])

# set up the points at which we want to know power spectrum and bispectrum in equilateral configuration
for ii in range(0,100):
    PhiExit = params[2] + .4 - ii*0.01
    PhiOut = np.append(PhiOut, PhiExit)
    k = MTS.kexitPhi(PhiExit, 1, back, params, MTSE)
    kOut= np.append(kOut, k)

zzOut , zzzOut, times = MTS.eqSpectra(kOut, back, params, 5.0,tols, MTSE)

fnlOut = 5./6*zzzOut/(3.0*zzOut**2)
############################################################################################################################################




################################################# Some plots ###############################################################################    
fig2 = plt.figure(2)
plt.plot(PhiOut, fnlOut, 'g',linewidth = 2)
title('Reduced bispectrum in equilateral configuration',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel('$fNL$', fontsize=20)
xlabel('$\phi^*$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.xlim(min(PhiOut),max(PhiOut))
plt.savefig("CEL1.png")


fig3 = plt.figure(3)
logk = np.log(kOut/kOut[-1])
plt.plot(logk, fnlOut, 'g',linewidth = 2)
title('Reduced bispectrum in equilateral configuration',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$fNL$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.xlim(np.min(logk),np.max(logk))
plt.savefig("CEL2.png")




fig4 = plt.figure(4)
G = 9./10*5./6.*1/3.*zzzOut*kOut**6/zzOut[-1]**2/kOut[-1]**6
plt.plot(logk, G, 'g',linewidth = 2)
title(r'G quantity',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$G/k^3$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.xlim(np.min(logk),np.max(logk))
plt.savefig("CEL3.png")


fig5 = plt.figure(5)
plt.plot(logk, np.log(zzOut/zzOut[-1]*kOut**3/kOut[-1]**3), 'g',linewidth = 2)
title('Power spectrum',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$\log({\cal P}/{\cal P}_{\rm pivot})$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.xlim(np.min(logk),np.max(logk))
plt.savefig("CEL4.png")


plt.show(fig1)
plt.show(fig2)
plt.show(fig3)
plt.show(fig4)
plt.show(fig5)