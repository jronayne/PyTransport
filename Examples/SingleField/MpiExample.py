#################### generate equilateral bispectrum using MTeasyPyStep setup ############################################################

from matplotlib import pyplot as plt

from pylab import *
import sys
import math 
import numpy as np

from mpi4py import MPI

location = "/Users/David/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransStep as PyT  # import module
import PyTransScripts as PyS

comm = MPI.COMM_WORLD
############################################################################################################################################




########################### initial field values ###########################################################################################

nF=PyT.nF() # gets number of fields (useful check)
nP=PyT.nP() # gets number of parameters needed (useful check)

fields = np.array([17.0])
tols = np.array([10**-8,10**-8])

params = np.zeros(nP)
params[0]=pow(10.,-5); params[1]=0.0018; params[2]=14.84; params[3]=0.022;

V = PyT.V(fields,params) # calculate potential from some initial conditions
dV=PyT.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)

initial = np.array([fields,-dV/np.sqrt(3.*V)]) # set initial conditions to be in slow roll
############################################################################################################################################




################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 40.0
t=np.linspace(Nstart, Nend, 1000)
back = PyT.backEvolve(t, initial, params,tols,False)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
############################################################################################################################################




################################## Calculate power spectrum and bispectrum for equilateral configs ########################################
points = 1000
#point = 100
kOutL=np.array([])
PhiOutL = np.array([])

rank=comm.Get_rank()
size=comm.Get_size()
num = points/size;

if rank==0 and float(points)/size != float(points/size):
    print "\n \n \n warning! number of points is divisable by number of processes, range will be reduced from that expected \n \n \n "


for ii in range(rank*num,rank*num+num):
    
    PhiExit = params[2] + .36 - ii*0.001
        
    
    #-0.3 + ii*0.011
    k = PyS.kexitPhi(PhiExit, 1, back, params, PyT)
    kOutL= np.append(kOutL, k)
    PhiOutL = np.append(PhiOutL, PhiExit)



zzL, zzzL, timesL = PyS.eqSpectra(kOutL, back, params, 5.0, tols, PyT)




########################################## Consolidate data on root process ################################################################
comm.Send(zzzL,dest=0)
comm.Send(zzL,dest=0)
comm.Send(timesL,dest=0)
comm.Send(kOutL,dest=0)
comm.Send(PhiOutL,dest=0)

if rank == 0:
    zzzOut = np.array([])    
    zzOut = np.array([])            
    timesOut = np.array([])    
    PhiOut = np.array([])    
    kOut   = np.array([])     
    for jj in range(0,size):    
        comm.Recv(zzzL,source = jj)
        comm.Recv(zzL,source = jj)                
        comm.Recv(timesL,source = jj)                
        comm.Recv(kOutL,source = jj)                
        comm.Recv(PhiOutL,source = jj)                
        
        zzzOut=np.append(zzzOut,zzzL)
        zzOut=np.append(zzOut,zzL)
        timesOut=np.append(timesOut,timesL)
        PhiOut=np.append(PhiOut,PhiOutL)
        kOut=np.append(kOut,kOutL)

    fnlOut = 5.0/6*zzzOut/(3.0*zzOut**2)
############################################################################################################################################





##############################################   Plots   ######################################################################################
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

    np.savetxt('test.dat', (kOut,PhiOut,fnlOut,G,zzOut,timesOut))   # x,y,z equal sized 1D arrays


    plt.show(fig2)
    plt.show(fig3)
    plt.show(fig4)
    plt.show(fig5)