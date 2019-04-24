#This file is part of PyTransport.
    
#PyTransport is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
    
#PyTransport is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
    
#You should have received a copy of the GNU General Public License
#along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.


# python code contains some useful scripts to use with the compiled PyTrans module.

import numpy as np
from scipy import interpolate 
import timeit
import sys

# this script finds initial conditions at least NBMassless e-folds before the massless point
# back must finely sample the backgroudn evolution for the initial conditions to be close to exactly NBMassless before

def unPackAlp(threePtOut,MTE):
    nF = MTE.nF()
    if np.size(threePtOut[0,:])!=1+4+2*nF+6*2*nF*2*nF+2*nF*2*nF*2*nF:
        print ("\n\n\n\n warning array you asked to unpack is not of correct dimension \n\n\n\n")
        return np.nan, np.nan, np.nan, np.nan, np.nan
    zetaMs = threePtOut[:,1:5]
    sig1R  = threePtOut[:,1+4+2*nF:1+4+2*nF+2*nF*2*nF]
    sig2R  = threePtOut[:,1+4+2*nF+2*nF*2*nF:1+4+2*nF+2*2*nF*2*nF]
    sig3R  = threePtOut[:,1+4+2*nF+2*2*nF*2*nF:1+4+2*nF+3*2*nF*2*nF]
    alp  = threePtOut[:,1+4+2*nF+6*2*nF*2*nF:]

    return zetaMs, np.reshape(sig1R,(np.size(threePtOut[:,0]), 2*nF,2*nF)), np.reshape(sig2R,(np.size(threePtOut[:,0]), 2*nF,2*nF)), np.reshape(sig3R,(np.size(threePtOut[:,0]), 2*nF,2*nF)), np.reshape(alpha,(np.size(threePtOut[:,0]), 2*nF,2*nF,2*nF))

def unPackSig(twoPtOut,MTE):
    nF = MTE.nF()
    if np.size(twoPtOut[0,:])!=1+1+2*nF+2*nF*2*nF:
        print ("\n\n\n\n warning array you asked to unpack is not of correct dimension \n\n\n\n")
        return np.nan , np.nan
    zeta = twoPtOut[:,1]
    sig  = twoPtOut[:,1+1+2*nF:1+1+2*nF+2*nF*2*nF]

    return zeta, np.reshape(sig,(np.size(twoPtOut[:,0]), 2*nF,2*nF))



def ICsBM(NBMassless,k, back, params, MTE):
    nF = np.size(back[0,1:])//2
    massEff = -1;    
    
    # calculate the element of back for which -k^2/a^2 + M^2 for M the largest eigenvalue of mass matrix
    jj=0
    while (massEff <0 and jj<np.size(back[:,0])-1):
        w, v = np.linalg.eig(MTE.ddV(back[jj,1:1+nF],params))
        eigen = np.max(w)
        massEff = -k**2*np.exp(-2.0*back[jj,0]) + eigen
        jj=jj+1
    if jj == np.size(back[:,0]):
        print ("\n\n\n\n warning massless condition not found \n\n\n\n")
        return np.nan, np.nan
    NMassless = back[jj-2,0]
    backExitMinus = np.zeros(2*nF)
    ll=0
    Ncond=-1.
    while (Ncond <0.0 and ll<np.size(back[:,0])-1):
        Ncond = back[ll,0]-(NMassless-NBMassless)
        ll=ll+1
        
    if ll == np.size(back[:,0]) or (NMassless-back[0,0])<NBMassless:
        print ("\n\n\n\n warning initial condition not found \n\n\n\n")
        return np.nan, np.nan

    NexitMinus = back[ll-2,0]
    backExitMinus = back[ll-2,1:]
        
    return NexitMinus, backExitMinus

# this script finds initial conditions at least NBExit e-folds before horizon exit of k
# back must finely sample the backgroudn evolution for the initial conditions to be close to exactly NBMassless before

def ICsBE(NBExit,k, back, params, MTE):
    nF = np.size(back[0,1:])//2
    kvaH = -1.;    
    jj=0
    while (kvaH <0.0 and jj<np.size(back[:,0])-1):
        H=MTE.H(back[jj,1:1+2*nF],params)
        kvaH = -k+np.exp(back[jj,0])*H
        jj=jj+1
    if jj == np.size(back[:,0]):
        print ("\n\n\n\n warning exit condition not found \n\n\n\n")
        return np.nan, np.nan
    NExit = back[jj-2,0]
    ll=0
    Ncond=-1
    while (Ncond <0 and ll<np.size(back[:,0])-1):
        Ncond = back[ll,0]-(NExit-NBExit)
        ll=ll+1
        
    if ll == np.size(back[:,0]) or (NExit-back[0,0])<NBExit:
        print  ("\n\n\n\n warning initial condition not found \n\n\n\n")
        return np.nan, np.nan
    NexitMinus = back[ll-2,0]
    backExitMinus = back[ll-2,1:]
        
    return NexitMinus, backExitMinus    
    

# find the earliest condition between the massless one and the horizon exit one
def ICs(NB, k, back, params, MTE):
    NBEs, fieldBE = ICsBE(NB,k,back, params, MTE)
    NBMs, fieldBM = ICsBM(NB,k,back, params, MTE)
    
    if (NBEs < NBMs) :
        return NBEs, fieldBE
    return NBMs, fieldBM


# calculates the power spectrum at each element in kA at the end of the background evolution (back)
def pSpectra(kA, back, params, NB, tols, MTE):
    zzOut=np.array([])
    times = np.array([])
    num = np.size(kA)
    
    for ii in range(0,num):
        print ("\n \n \n performing " + str(ii+1) + " of " + str(num)+"\n \n \n")
        k=kA[ii]
        Nstart, backExitMinus = ICs(NB, k, back, params, MTE)
        start_time = timeit.default_timer()

        if Nstart == np.nan:
            twoPt = numpy.empty((2,2))
            twoPt[:] = np.nan
        else:
            t=np.linspace(Nstart,back[-1,0], 10)
            # run solver for this triangle
            twoPt = MTE.sigEvolve(t,k, backExitMinus,params,tols,True) # all data from three point run goes into threePt array
        zzOut=np.append(zzOut, twoPt[-1,1])
        times = np.append(times, timeit.default_timer()-start_time)
        
    return zzOut, times

# calculates the power spectrum at each element in kA at the end of the background evolution (back) in a manner suitable to be called over many processes
def pSpecMpi(kA, back, params, NB, tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    rank=comm.Get_rank()
    size=comm.Get_size()
    points  = np.size(kA)
    num = points/size;

    if  float(points)/size != float(points//size):
        if rank==0:
            print ("\n \n \n warning! number of points is divisable by number of processes, exiting \n \n \n ")
        return (np.empty, np.empty)

    kOutL = kA[rank*num:rank*num+num]
    
    zzL, timesL = pSpectra(kOutL, back, params, NB, tols, MTE)

    if rank != 0:
        comm.Send(zzL,dest=0)
        comm.Send(timesL,dest=0)

    if rank == 0:
        zzzOut = np.array([])    
        zzOut = np.array([])            
        timesOut = np.array([])      
        zzOut=np.append(zzOut,zzL)
        timesOut=np.append(timesOut,timesL)

        for jj in range(1,size):
            comm.Recv(zzL,source = jj)                
            comm.Recv(timesL,source = jj)                
            zzOut=np.append(zzOut,zzL)
            timesOut=np.append(timesOut,timesL)

        return (zzOut, timesOut)
    else:
        return (np.empty, np.empty)

# calculates the power spectrum and bisecpturm in equilateral configuration at each element in kA at the end of the background evolution (back)
def eqSpectra(kA, back, params, NB, tols, MTE):
    zzzOut=np.array([])
    zzOut=np.array([])
    times = np.array([])
    num = np.size(kA)
    
    for ii in range(0,num):
        print ("\n \n \n performing " + str(ii+1) + " of " + str(num) + "\n \n \n")
        k=kA[ii]
        Nstart, backExitMinus = ICs(NB, k, back, params, MTE)
        t=np.linspace(Nstart,back[-1,0], 10)
        k1=k;k2=k;k3=k;
        # run solver for this triangle
        start_time = timeit.default_timer()
        
        if Nstart == np.nan:
            nF = MTE.nF(); threePt = numpy.empty((2,5))
            threePt[:] = np.nan
        else:
            threePt = MTE.alphaEvolve(t,k1,k2,k3, backExitMinus,params,tols, True) # all data from three point run goes into threePt array
        zzOut=np.append(zzOut, threePt[-1,1])
        zzzOut=np.append(zzzOut, threePt[-1,4])
        times = np.append(times, timeit.default_timer()-start_time)
    return zzOut, zzzOut, times

# calculates the power spectrum and bisecpturm in equilateral configuration at each element in kA at the end of the background evolution (back) in a manner suitable to be run accoss many processes
def eqSpecMpi(kA, back, params, NB,tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD

    rank=comm.Get_rank()
    size=comm.Get_size()
    points  = np.size(kA)
    num = points/size;

    if float(points)/size != float(points//size):
        if rank==0:
            print ("\n \n \n warning! number of points is divisable by number of processes, exiting \n \n \n ")
        return (np.empty, np.empty, np.empty)

    kOutL = kA[rank*num:rank*num+num]

    zzL, zzzL, timesL = eqSpectra(kOutL, back, params, NB, tols, MTE)
  
    if rank !=0:
        comm.Send(zzzL,dest=0)
        comm.Send(zzL,dest=0)
        comm.Send(timesL,dest=0)

    if rank == 0:
        zzzOut = np.array([])    
        zzOut = np.array([])            
        timesOut = np.array([])      
        zzzOut=np.append(zzzOut,zzzL)
        zzOut=np.append(zzOut,zzL)
        timesOut=np.append(timesOut,timesL)
        
        for jj in range(1,size):
            comm.Recv(zzzL,source = jj)
            comm.Recv(zzL,source = jj)                
            comm.Recv(timesL,source = jj)                
        
            zzzOut=np.append(zzzOut,zzzL)
            zzOut=np.append(zzOut,zzL)
            timesOut=np.append(timesOut,timesL)
        
        return (zzOut, zzzOut, timesOut)
    else: 
        return (np.empty, np.empty, np.empty)

# calcualtes the bispectrum in the alpha beta notation for a given kt at every value of the alphaIn and betIn arrays. The bispectrum is given an nsnaps times always incuding the final time of the evolution (back)
def alpBetSpectra(kt,alphaIn, betaIn, back, params, NB, nsnaps, tols, MTE):

    Hin = np.zeros(np.size(back[:,0]))
    for jj in range(0,np.size(back[:,0])):
            Hin[jj]=MTE.H(back[jj,1:],params)
    aH=np.exp(back[:,0])*Hin
    positions = np.argsort(aH)
    Nexit = interpolate.splev(kt/3., interpolate.splrep(aH[positions], back[positions,0], s=1e-15),der=0)
    Nend = back[-1,0]
    snaps =  np.linspace(Nexit-(NB - .1),Nend,nsnaps)
    if (nsnaps ==1 or nsnaps == 0):
        snaps= np.array([Nend]) 
        nsnaps =1
          
    biAOut = np.zeros([np.size(alphaIn),np.size(betaIn),np.size(snaps)])
    zz1 = np.zeros([np.size(alphaIn),np.size(betaIn),np.size(snaps)])
    zz2 = np.zeros([np.size(alphaIn),np.size(betaIn),np.size(snaps)])
    zz3 = np.zeros([np.size(alphaIn),np.size(betaIn),np.size(snaps)])
    times = np.zeros([np.size(alphaIn),np.size(betaIn)])
    
    for l in range(0,np.size(alphaIn)):
      alpha =  alphaIn[l]          
      for j in range(0,np.size(betaIn)):
         print ("\n \n \n performing " + str(l+1)+ " " +str(j+1) + " of " + str(np.size(alphaIn))+ " " + str(np.size(betaIn)) + "\n \n \n")
         timebefore =     timeit.default_timer()
         beta =betaIn[j]
         k1 = kt/2. - beta*kt/2.
         k2 = kt/4.*(1.+alpha+beta)
         k3 = kt/4.*(1.-alpha+beta)
         if alpha > -(1-beta) and alpha < 1-beta:
            kM= min(k1,k2,k3)
            Nstart, backExitMinus = ICs(NB, kM, back, params, MTE)
        # run solver for this triangle
            t= np.concatenate((np.array([Nstart]),snaps))
            if Nstart == np.nan:
                    threePt = numpy.empty((2,5))
                    threePt[:] = np.nan
            else:
                    threePt = MTE.alphaEvolve(t,k1,k2,k3, backExitMinus,params,tols,True)
            zzz= threePt[:,:5]
            
            for ii in range(1,nsnaps+1):
                biAOut[l, j, ii-1]=zzz[ii,4]
                zz1[l, j, ii-1] = zzz[ii,1]
                zz2[l, j, ii-1] = zzz[ii,2]
                zz3[l, j, ii-1] = zzz[ii,3]
         else:
            for ii in range(0,nsnaps+1):
                biAOut[l, j, ii-1]=np.nan
                zz1[l, j, ii-1] = np.nan
                zz2[l, j, ii-1] = np.nan
                zz3[l, j, ii-1] = np.nan
         
         times[l,j] = timeit.default_timer()-timebefore 
    return (biAOut, zz1, zz2, zz3, times, snaps)

# performs the same task as alpBetSpectra but in a manner suitable to be spread across many processes
def alpBetSpecMpi(kt, alpha, beta, back, params, NB, nsnaps,tols, MTE):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    side = np.size(alpha)
    Nbefore=NB
    rank=comm.Get_rank()
    size=comm.Get_size()
    num = side/size;
    
    if float(side)/size != float(side//size):
        if rank ==0:
            print ("\n \n \n warning! number size of alpha must be divisable by number of processes, exiting \n \n \n ")
        return (np.empty, np.empty, np.empty, np.empty, np.empty)
    else:

        alphaL=alpha[rank*num:rank*num+num]
    
        BzL, Pz1L, Pz2L, Pz3L,  timesL, snaps = alpBetSpectra(kt,alphaL, beta, back, params, Nbefore, nsnaps,tols, MTE)

        if rank != 0:
            comm.send(Pz1L,dest=0)
            comm.send(Pz2L,dest=0)
            comm.send(Pz3L,dest=0)
            comm.send(BzL,dest=0)
            comm.send(timesL,dest=0)

        if rank == 0:
            Bztot = np.zeros([np.size(alpha),np.size(beta),np.size(BzL[0,0,:])])
            Pz1tot = np.zeros([np.size(alpha),np.size(beta),np.size(BzL[0,0,:])])
            Pz2tot = np.zeros([np.size(alpha),np.size(beta),np.size(BzL[0,0,:])])
            Pz3tot = np.zeros([np.size(alpha),np.size(beta),np.size(BzL[0,0,:])])
            timestot = np.zeros([np.size(alpha),np.size(beta)])
            Bztot[0:num,:,:]=BzL
            Pz1tot[0:num,:,:]=Pz1L
            Pz2tot[0:num,:,:]=Pz2L
            Pz3tot[0:num,:,:]=Pz3L
            timestot[0:num,:]=timesL
            
            for jj in range(1,size):
                Pz1L=comm.recv(source=jj)
                Pz2L=comm.recv(source=jj)
                Pz3L=comm.recv(source=jj)
                BzL=comm.recv(source = jj)
                timesL = comm.recv(source = jj)
        
                Bztot[jj*num:jj*num+num,:,:]=BzL
                Pz1tot[jj*num:jj*num+num,:,:]=Pz1L
                Pz2tot[jj*num:jj*num+num,:,:]=Pz2L
                Pz3tot[jj*num:jj*num+num,:,:]=Pz3L
                timestot[jj*num:jj*num+num,:]=timesL
            return (Bztot, Pz1tot, Pz2tot, Pz3tot,  timestot, snaps)
        else:
            return (np.empty, np.empty, np.empty, np.empty, np.empty, snaps)



def kexitN(Nexit, back, params, MTE, exact=False):
    nF = np.size(back[0,1:])//2  
    backExit = np.zeros(2*nF)

    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
    k = np.exp(Nexit) *  MTE.H(backExit,params);
    if exact==True:
        indx = np.argmin(np.abs(np.array(back[:,0])-Nexit))#
        k = np.exp(Nexit) *  MTE.H(back[indx,1:],params);

    return k


def kexitPhi(PhiExit, n, back, params, MTE):
    nF = np.size(back[0,1:])//2    
    backExit = np.zeros(2*nF)
    positions = np.argsort(back[:,n])
    Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[positions,n], back[positions,0], s=1e-15),der=0)
    
    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
    k = np.exp(Nexit) *  MTE.H(backExit,params);
    return k
    

