#################### generate alpha beta  bispectrum using PseudoScalar and MPI ############################################################
from matplotlib import pyplot as plt   # import package for plotting
from mpl_toolkits.mplot3d import Axes3D
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below

from mpi4py import MPI

location = "/home/jwr/Code/June/PyTransport2Dist/PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTransport uses

import PyTransPseudo as PyT  # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD
tols = np.array([10**-8,10**-8])

########################### initial field values ###########################################################################################
fields = np.array([10.0,0.01,13.0]) # we set up a numpy array which contains the values of the fields

nP=PyT.nP()   # the .nP() function gets the number of parameters needed for the potential -- this can be used as a useful crosscheck
pvalue = np.zeros(nP)
M=1.*pow(10.,-6.0)

pvalue[0]=1./2. *30.0 * M**2; pvalue[1]=1./2. *300.0 * M**2; pvalue[2]=1./2. *30./81. *M**2; # we set up numpy array which contains values of the parameters

nF=PyT.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)

V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV=PyT.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)

initial = np.concatenate((fields,np.array([0,0,0]))) # sets an array containing field values and there derivative in cosmic time
                                                      # (set using the slow roll equation)
############################################################################################################################################


################################## run the background fiducial run ########################################################################
Nstart = 0.0
Nend = 70.0
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, pvalue,tols,False)
###########################################################################################################################################

side = 100
nsnaps = 0
Nbefore=4.5
rank=comm.Get_rank()

NExit = 60.0
kt = PyS.kexitN(NExit, back, pvalue, PyT)

alpha =  np.linspace(-1,1,side)
beta=np.linspace(0,1,side/2)

Bztot, Pz1tot, Pz2tot, Pz3tot,  times, snaps = PyS.alpBetSpecMpi(kt,alpha, beta, back, pvalue, Nbefore, nsnaps,tols, PyT)

if rank == 0:
    bet, alp = np.meshgrid(beta, alpha)
    np.save('data15/alp',alp);np.save('data15/bet',bet); np.save('data15/alBetBi',Bztot)
    np.save('data15/alBetPz1',Pz1tot); np.save('data15/alBetPz2.npy',Pz2tot); np.save('data15/alBetPz3',Pz3tot)
print ("\n\n process " + str(rank) + " done \n\n")