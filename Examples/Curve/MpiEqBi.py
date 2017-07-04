#################### generate equilateral bispectrum using the quasi-single field and MPI ############################################################
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

import PyTransCurve as PyT  # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD
tols = np.array([10**-10,10**-10])

########################### set some field values and field derivatives in cosmic time ####################################################
omega = pi/30.0
R0 = np.sqrt(10.0**(-10)/3.0) / (omega *np.sqrt(10.0**(-9)))

fields = np.array([-R0, (1e-2)*R0]) # we set up a numpy array which contains the values of the fields

nP=PyT.nP()   # the .np function gets the number of parameters needed for the potential -- this can be used as a useful crosscheck
pvalue = np.zeros(nP)
pvalue[1]=1./np.sqrt(3.0); pvalue[2]=1.0/10.0**(-5); pvalue[3]=1.0/omega**.5 * 1./2. * (10.0**(-10))**(-3./4.)
pvalue[0]=R0
# we set up numpy array which contains values of the parameters

nF=PyT.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)

V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV = PyT.dV(fields,pvalue) # calculate derivatives of potential

initial = np.concatenate((fields, np.array([0,0]))) # sets an array containing field values and there derivative in cosmic time
# (set using the slow roll equation)

############################################################################################################################################



################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 28.0
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = PyT.backEvolve(t, initial, pvalue,tols,False) # The output is read into the back numpy array

############################################################################################################################################

rank=comm.Get_rank()
points = 500
Nexit1=17.0
kOut =np.array([])
NOut = np.array([])
for ii in range(0,points):
    Nexit = Nexit1+0.014*ii
    k = PyS.kexitN(Nexit, back, pvalue, PyT)
    kOut= np.append(kOut,k)
    NOut = np.append(NOut,Nexit)

Pztot, Bztot, times = PyS.eqSpecMpi(kOut, back, pvalue, 5.0, tols,PyT)

print ("\n\n process " + str(rank) + " done \n\n")

if rank ==0:
    fnlOut = 5.0/6*Bztot/(3.0*Pztot**2.0)
    plt.plot(NOut, fnlOut, linewidth=2)
    #plt.plot(np.log(kOut/kOut[0]), fnlOut, linewidth=2)
    title('Reduced bispectrum in equilateral configuration',fontsize=15) ;grid(True); plt.legend(fontsize=15); ylabel('$fNL$', fontsize=20)
    xlabel(r'e-fold exit time', fontsize=15); grid(True); plt.legend(fontsize=15); 
    #xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15); grid(True); plt.legend(fontsize=15); 
    #plt.xlim(min(np.log(kOut/kOut[0])),max(np.log(kOut/kOut[0])));
    plt.savefig("BiEq.png")
    
    np.savetxt('data/EqBi.dat', (kOut,NOut,fnlOut,Bztot,Pztot))   # x,y,z equal sized 1D arrays
