#################### generate equilateral bispectrum using PyTransDQuad and MPI ############################################################
from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below

from mpi4py import MPI

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransDQuad as PyT  # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD

########################### initial field values ###########################################################################################
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

tols = np.array([10**-8,10**-8])

################################## run the background fiducial run ########################################################################
Nstart = 0.0
Nend = 60.0
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, pvalue,tols,False)
###########################################################################################################################################

rank=comm.Get_rank()
points = 500

NExit1 = 10.0
NExit2 = 10.0+0.08*points
k1 = PyS.kexitN(NExit1, back, pvalue, PyT)
k2 = PyS.kexitN(NExit2, back, pvalue, PyT)
kOut= np.logspace(log10(k1), log10(k2), points)

Pztot, Bztot, times = PyS.eqSpecMpi(kOut, back, pvalue, 4.5, tols,PyT)

print ("\n\n process " + str(rank) + " done \n\n")

if rank ==0:
    fnlOut = 5.0/6*Bztot/(3.0*Pztot**2.0)
    fig2 = plt.figure(2)
    plt.plot(np.log(kOut/kOut[0]), fnlOut, linewidth=2)
    title('Reduced bispectrum in equilateral configuration',fontsize=15) ;grid(True); plt.legend(fontsize=15); ylabel('$fNL$', fontsize=20)
    xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15); grid(True); plt.legend(fontsize=15); plt.xlim(min(np.log(kOut/kOut[0])),max(np.log(kOut/kOut[0]))); 
    plt.savefig("BiEq.png")
    
    np.savetxt('data/EqBi.dat', (kOut,fnlOut,Bztot,Pztot))   # x,y,z equal sized 1D arrays
