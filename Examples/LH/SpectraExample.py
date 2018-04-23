####################################### Langlois Spectra example of basic PyTransport functions ###########################################
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

#This file contains simple examples of using the PyTransport package for the heavy field example of Langlois.
#It assumes the LangHeavySeptup file has been run to install a LH version of PyTransport
#It is recommended you restart the kernel to insure any updates to PyTransLH are imported 

location = "/Users/mulryne/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # his add sets the other paths that PyTransport uses

import PyTransLH as PyT;  # import module  
import PyTransScripts as PyS;

# Example 

########################### set initial field values and parameters for a simple example run ###################################################
shift=231.
fields = np.array([-2.0-100.*math.sqrt(6.) +shift, 2.0*math.tan(math.pi/20.)])

nF=PyT.nF() # gets number of fields (useful check)
nP=PyT.nP() # gets number of parameters needed (useful check)

params = np.zeros(nP)
params[0]=1.*pow(10.,-7); params[1]=1.*pow(10.,-4); params[2]=math.pi/10.; params[3] = -100.*math.sqrt(6.) +shift; params[4]=10.*math.sqrt(3.);

V = PyT.V(fields,params) # calculate potential from some initial conditions
dV=PyT.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions using slow roll expression  
############################################################################################################################################
###########################################################################################################################################
tols = np.array([10**-8,10**-8])

################################## run the background fiducial run #########################################################################
Nstart = -5.0
Nend = 30.
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, params,tols,False)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g') # always good to inspect this plot to make sure its senisble
############################################################################################################################################


fnlOut=np.array([])
kOut=np.array([])
PhiOut=np.array([])

# set up the points at which we want to know power spectrum and bispectrum in equilateral configuration
for ii in range(0,50):
    PhiExit = -100.*math.sqrt(6.) + shift  -0.1 + ii*0.011
    PhiOut = np.append(PhiOut, PhiExit)
    k = PyS.kexitPhi(PhiExit, 1, back, params, PyT) 
    kOut= np.append(kOut, k)

zzOut , zzzOut, times = PyS.eqSpectra(kOut, back, params, 4.0,tols, PyT)

fnlOut = 5./6*zzzOut/(3.0*zzOut**2)

    
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
plt.savefig("LH1.png")


fig3 = plt.figure(3)
plt.plot(np.log(kOut/kOut[-1]), fnlOut, 'g',linewidth = 2)
title('Reduced bispectrum in equilateral configuration',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$fNL$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
#plt.xlim((0,5.65))
plt.savefig("LH2.png")



fig4 = plt.figure(4)
G = 9./10*5./6.*1/3.*zzzOut*kOut**6/zzOut[-1]**2/kOut[-1]**6
plt.plot(log(kOut/kOut[-1]), G, 'g',linewidth = 2)
title(r'G quantity',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$G/k^3$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
#plt.xlim((0,5.65))
plt.savefig("LH3.png")


fig5 = plt.figure(5)
plt.plot(np.log(kOut/kOut[-1]), np.log(zzOut/zzOut[-1]*kOut**3/kOut[-1]**3), 'g',linewidth = 2)
title('Power spectrum',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$\log({\cal P}/{\cal P}_{\rm pivot})$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
#plt.xlim((0,5.65))
plt.savefig("LH4.png")


plt.show(fig1)
plt.show(fig2)
plt.show(fig3)
plt.show(fig4)
plt.show(fig5)
