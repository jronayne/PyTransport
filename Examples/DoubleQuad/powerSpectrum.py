#################### generate power spectrum using PyTransDQuad installed using the setup file which accompanies this one ####################
from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below
############################################################################################################################################


############################################################################################################################################
#This file contains simple examples of using the PyTransDQuad 
#It assumes the DQuadSetup file has been run to install a double quadratic version of PyTrans
#It is recommended you restart the kernel before running this file to insure any updates to PyTransDQuad are imported

location = "/Users/David/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransDQuad as PyT  # import module as PyT (PyTransDQuad is quite long to type each time and it saves time to use a shorter name
                             # using a generic name PyT means the file can be more easily reused for a different example (once field values 
                             # etc are altered)
import PyTransScripts as PyS  # import the scripts module as PyS for convenience
###########################################################################################################################################


########################### set up initial conditions for background run ##################################################################
fields = np.array([12.0, 12.0]) # we set up a numpy array which contains the values of the fields

nP=PyT.nP()   # the .nP() function gets the number of parameters needed for the potential -- this can be used as a useful crosscheck
pvalue = np.zeros(nP) 
pvalue[0]=10.0**(-5.0); pvalue[1]=9.0*10.0**(-5) # we set up numpy array which contains values of the parameters

nF=PyT.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)

V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV=PyT.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)

initial = np.concatenate((fields, -dV/np.sqrt(3* V))) # sets an array containing field values and there derivative in cosmic time 
                                                      # (set using the slow roll equation)
############################################################################################################################################


################################## run the background fiducial run ################################################################
Nstart = 0.0
Nend = 70.0
t=np.linspace(Nstart, Nend, 1000)
tols = np.array([10**-8,10**-8])

back = PyT.backEvolve(t, initial, pvalue,tols,False)
###########################################################################################################################################


###########################################################################################################################################
# set up an array of k values to calculate power spectrum as a function of ks
kOut = np.array([])
for ii in range(0,500):
    NExit = 10.0 + ii*0.08
    k = PyS.kexitN(NExit, back, pvalue, PyT)  
    kOut = np.append(kOut,k) # this builds an array of ks associated with different NExit times from 10 to 50
Pz, times = PyS.pSpectra(kOut,back,pvalue,4.0,PyT) # this cacalcute P_z for this range of ks, using 5.0 e-folds of subhorizon evolution
#==============================================================================
# zz, zzz, timesB = PyS.eqSpectra(kOut, back, pvalue, 4.0, PyT)
#==============================================================================


fig1=plt.figure(1)
plt.plot(np.log(kOut/kOut[0]), Pz/Pz[0] *kOut**3/kOut[0]**3, linewidth = 2 )
title('Power spectrum',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'$\log({\cal P}/{\cal P}_{\rm pivot})$', fontsize=20); 
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15);yscale('log'); plt.xlim(min(np.log(kOut/kOut[0])),max(np.log(kOut/kOut[0])));plt.savefig("Pz.png")

#fnl=5.0/6.0*zzz/(3.0*zz**2.0)  
#plt.plot(np.log(kOut/kOut[0]), fnl )
##fig3=plt.figure(3)
#title(r'$f_{NL}$',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'$f_{NL}$', fontsize=20); 
#xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15); plt.xlim(min(np.log(kOut/kOut[0])),max(np.log(kOut/kOut[0]))); plt.savefig("fnl.png") 
 
############################################################################################################################################


############################################################################################################################################
# finally lets find n_s
from scipy.interpolate import UnivariateSpline
derivativeSp = UnivariateSpline(np.log(kOut/kOut[0]), np.log(Pz),k=4, s=1e-15).derivative()
ns=derivativeSp(np.log(kOut/kOut[0]))+4.0

fig2=plt.figure(2)
plt.plot(np.log(kOut/kOut[0]), ns , linewidth = 2)
title('Spectral index',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'$n_s$', fontsize=20); 
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15); plt.xlim(min(np.log(kOut/kOut[0])),max(np.log(kOut/kOut[0]))); plt.savefig("ns.png")
