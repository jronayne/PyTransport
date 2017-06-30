#################### generate alpha beta  bispectrum using PyTransDQuad and MPI ############################################################
from matplotlib import pyplot as plt   # import package for plotting
from mpl_toolkits.mplot3d import Axes3D
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below

from mpi4py import MPI

location = "/Users/mulryne/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransDQuad as PyT  # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD
tols = np.array([10**-8,10**-8])

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


################################## run the background fiducial run ########################################################################
Nstart = 0.0
Nend = 50.0
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, pvalue,tols,False)
###########################################################################################################################################

side = 100
nsnaps = 0
Nbefore=4.5
rank=comm.Get_rank()

NExit = 15.0
kt = PyS.kexitN(NExit, back, pvalue, PyT)

alpha =  np.linspace(-1,1,side)
beta=np.linspace(0,1,side/2)

Bztot, Pz1tot, Pz2tot, Pz3tot,  times, snaps = PyS.alpBetSpecMpi(kt,alpha, beta, back, pvalue, Nbefore, nsnaps, PyT)

if rank == 0:
    bet, alp = np.meshgrid(beta, alpha)
    np.save('data/alp',alp);np.save('data/bet',bet); np.save('data/alBetBi',Bztot)
    np.save('data/alBetPz1',Pz1tot); np.save('data/alBetPz2.npy',Pz2tot); np.save('data/alBetPz3',Pz3tot)
print ("\n\n process " + str(rank) + " done \n\n")
