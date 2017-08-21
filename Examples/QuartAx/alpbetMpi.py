#################### generate alpha beta  bispectrum using PyTransAxQrt ############################################################

from matplotlib import pyplot as plt

from pylab import *
import sys
import math 
import numpy as np

from mpi4py import MPI

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransAxQrt as PyT  # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD



########################### set initial field values and parameters for a simple example run ###################################################
nF=PyT.nF() # gets number of fields (useful check)
nP=PyT.nP() # gets number of parameters needed (useful check)

fields = np.array([23.5,.5-0.001])

params = np.zeros(nP)
params[0]=1.*pow(10.,-10); params[1]=1.; params[2]=25.0**2.0*params[0]/4.0/math.pi**2;

V = PyT.V(fields,params) # calculate potential from some initial conditions
dV=PyT.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions using slow roll expression
############################################################################################################################################
tols = np.array([10**-8,10**-8])


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 70.0
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = PyT.backEvolve(t, initial, params,tols,False) # The output is read into the back numpy array
############################################################################################################################################

rank=comm.Get_rank()

side = 140
nsnaps = 150
Nbefore=4.5

NExit = 14.0
kt = PyS.kexitN(NExit, back, params, PyT)
kt =3.*kt
alpha =  np.linspace(-1,1,side)
beta=np.linspace(0,1,side/2)

Bztot, Pz1tot, Pz2tot, Pz3tot,  times, snaps = PyS.alpBetSpecMpi(kt,alpha, beta, back, params, Nbefore, nsnaps,tols, PyT)

if rank == 0:
    bet, alp = np.meshgrid(beta, alpha)
    np.save('data/alp',alp);np.save('data/bet',bet); np.save('data/alBetBi',Bztot); np.save('data/times',times)
    np.save('data/alBetPz1',Pz1tot); np.save('data/alBetPz2.npy',Pz2tot); np.save('data/alBetPz3',Pz3tot)
    np.save('data/snaps',snaps)
print "\n\n process", rank, "done \n\n"
