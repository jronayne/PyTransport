#################### generate equilateral bispectrum using Langlois model and MPI ############################################################
from matplotlib import pyplot as plt

from pylab import *
import sys
import math 
import numpy as np

from mpi4py import MPI

location = "/Users/mulryne/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransLH as PyT   # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD

########################### initial field values ###########################################################################################
nP=PyT.nP()
shift=231.
fields = np.array([-2.-100.*math.sqrt(6.) +shift, 2.*math.tan(math.pi/20.)])
pvalue = np.zeros(nP)
pvalue[0]=1.*pow(10.,-7); pvalue[1]=1.*pow(10.,-4); pvalue[2]=math.pi/10.; pvalue[3] = -100.*math.sqrt(6.) +shift; pvalue[4]=10.*math.sqrt(3.);# gelaton-like case (slow turn)
#pvalue[0]=1.*pow(10.,-7); pvalue[1]=1.*pow(10.,-4); pvalue[2]=math.pi/10.; pvalue[3] = -100.*math.sqrt(6.) +shift; pvalue[4]=1000.*math.sqrt(3.);# oscillating  case (particle production case in paper)


nF=PyT.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV=PyT.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions to be in slow roll
############################################################################################################################################
tols = np.array([10**-8,10**-8])


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 40.
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, pvalue,tols,False)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
############################################################################################################################################

rank=comm.Get_rank()

points = 500

for ii in range(0,points):
    PhiExit = -100.*math.sqrt(6.) +shift -0.3 + ii*0.0015
    #PhiExit = -100.*math.sqrt(6.) +shift -0.3 + ii*0.0025 +1.0
    k = PyS.kexitPhi(PhiExit, 1, back, pvalue, PyT)
    kOut= np.append(kOut, k)
    PhiOut = np.append(PhiOut, PhiExit)

zztot, zzztot, timestot = PyS.eqSpecMpi(kOut, back, pvalue, 5.0,tols, PyT)

if rank ==0:
    fnlOut = 5.0/6*zzztot/(3.0*zztot**2)
    np.savetxt('dataT/EqBi.dat', (kOut,PhiOut,fnlOut,zzztot,zztot))   # x,y,z equal sized 1D arrays


    
    fig2 = plt.figure(2)
    plt.plot(phiOut, fnlOut, 'g',linewidth = 2)
    title('Reduced bispectrum in equilateral configuration',fontsize=15)
    grid(True)
    plt.legend(fontsize=15)
    ylabel('$fNL$', fontsize=20)
    xlabel('$\phi^*$', fontsize=15)
    grid(True)
    plt.legend(fontsize=15)
    plt.xlim(min(phiOut),max(phiOut))
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
    G = 9./10*5./6.*1/3.*zzztot*kOut**6/zztot[-1]**2/kOut[-1]**6
    plt.plot(np.log(kOut/kOut[-1]), G, 'g',linewidth = 2)
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
    plt.plot(np.log(kOut/kOut[-1]), np.log(zztot/zztot[-1]*kOut**3/kOut[-1]**3), 'g',linewidth = 2)
    title('Power spectrum',fontsize=15)
    grid(True)
    plt.legend(fontsize=15)
    ylabel(r'$\log({\cal P}/{\cal P}_{\rm pivot})$', fontsize=20)
    xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
    grid(True)
    plt.legend(fontsize=15)
    #plt.xlim((0,5.65))
    plt.savefig("LH4.png")
    
    np.savetxt('test.dat', (kOut,phiOut,fnlOut,G,zztot))   # x,y,z equal sized 1D arrays
    
    

    plt.show(fig2)
    plt.show(fig3)
    plt.show(fig4)
    plt.show(fig5)