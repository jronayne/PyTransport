#################### generate equilateral bispectrum using PyTransStep and MPI ############################################################
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below

from mpi4py import MPI

location = "/Users/mulryne/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransStep as PyT  # import module
import PyTransScripts as PyS

########################### initial field values ###########################################################################################

nF=PyT.nF() # gets number of fields (useful check)
nP=PyT.nP() # gets number of parameters needed (useful check)

fields = np.array([16.5])

pvalue = np.zeros(nP)
pvalue[0]=pow(10.,-5); pvalue[1]=0.0018; pvalue[2]=14.84; pvalue[3]=0.022;

V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV=PyT.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)

initial = np.array([fields,-dV/np.sqrt(3.*V)]) # set initial conditions to be in slow roll
############################################################################################################################################


################################## run the background fiducial run #########################################################################
tols = np.array([10**-8,10**-8])

Nstart = 0.0
Nend = 40.0
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, pvalue,tols,False)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
############################################################################################################################################

side = 140
nsnaps = 0
Nbefore=4.5

#PhiExit = 14.6
Nexit =14.8
k = PyS.kexitN(Nexit, back, pvalue, PyT)
kt = 3*k

alpha =  np.linspace(-1,1,side)
beta=np.linspace(0,1,side/2)

Bztot, Pz1tot, Pz2tot, Pz3tot,  times = PyS.alpBetSpecMpi(kt,alpha, beta, back, pvalue, Nbefore, nsnaps,tols, PyT)

comm = MPI.COMM_WORLD
rank =comm.Get_rank()
if rank == 0:

    fnlOut = 5.0/6*Bztot[:,:,-1]/(Pz1tot[:,:,-1]*Pz2tot[:,:,-1] + Pz1tot[:,:,-1]*Pz3tot[:,:,-1] + Pz2tot[:,:,-1]*Pz3tot[:,:,-1])
    bet, alp = np.meshgrid(beta, alpha)

    np.save('data/alp',alp);np.save('data/bet',bet); np.save('data/alBetBi',Bztot)
    np.save('data/alBetPz1',Pz1tot); np.save('data/alBetPz2.npy',Pz2tot); np.save('data/alBetPz3',Pz3tot)
print "\n\n process", rank, "done \n\n"
