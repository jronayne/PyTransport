####################################### PyTransPyStep simple example of basic functions ###########################################

from matplotlib import pyplot as plt
import time
import math 
import numpy as np
import sys 




############################################################################################################################################

#This file contains simple examples of using the PyTrans package for the single field Step example of Chen et al.
#It assumes the StepPotSetup file has been run to install a Step version of PyTrans
#It is recommended you restart the kernel before running this file to insure any updates to PyTransPyStep are imported

location = "/Users/David/Dropbox/MTeasyDist/MTeasy/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # his add sets the other paths that MTeasy uses

import PyTransPyStep as MTSE  # import module
import PyTransScripts as MTS
###########################################################################################################################################

# Example 

########################### set initial field values and parameters for a simple example run ###############################################
tols = np.array([10**-8,10**-8])

nF=MTSE.nF() # gets number of fields (useful check)
nP=MTSE.nP() # gets number of parameters needed (useful check)


params = np.zeros(nP)
params[0]=pow(10.,-5); params[1]=0.0018; params[2]=14.84; params[3]=0.022;

fields = np.array([17.0])

V = MTSE.V(fields,params) # calculate potential from some initial conditions
dV=MTSE.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)

initial = np.array([fields,-dV/np.sqrt(3.*V)]) # set initial conditions to be in slow roll
############################################################################################################################################






################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 40.0
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = MTSE.backEvolve(t, initial, params,tols,False) # The output is read into the back numpy array

# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,1 ], 'g')
############################################################################################################################################






############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t

PhiExit = params[2]
k = MTS.kexitPhi(PhiExit, 1, back, params, MTSE) 

# other scales can then be defined wrt to k
############################################################################################################################################






################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds of massless sub-horizon evolution using a spline routine
NBMassless = 4.0
Nstart, backExitMinus = MTS.ICsBM(NBMassless, k, back, params, MTSE)
print "NStart= ", Nstart 
print k
t=np.linspace(Nstart,Nend, 1000)  # array at which output is returned -- itial valye should correspond to initial field values
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an crude estimate for ns

twoPt = MTSE.sigEvolve(t, k, backExitMinus,params, tols,True) # puts information about the two point fuction in twoPt array
zz1=twoPt[:,1] # the second column is the 2pt of zeta
zz1a=zz1[-1] 
twoPt=MTSE.sigEvolve(t, k+.1*k, backExitMinus,params, tols,True)
zz2=twoPt[:,1]
zz2a=zz2[-1]
n_s = (np.log(zz2a)-np.log(zz1a))/(np.log(k+.1*k)-np.log(k))+4.0
print n_s
fig2 = plt.figure(2)
plt.plot(t, zz1,'r')
###########################################################################################################################################






###################################### example bispectrum run ##############################################################################
# set three scales in FLS manner (using alpha beta notation)
alpha=0.
beta =1/3.

k1 = k/2 - beta*k/2.
k2 = k/4*(1+alpha+beta)
k3 = k/4*(1-alpha+beta)
k1 =k; k2 =k ;k3 =k;
# find initial conditions for 4 e-folds of massless evolution for the largest k (which stays inside the longest)

kM = np.max(np.array([k1,k2,k3]))
Nstart, backExitMinus = MTS.ICsBM(NBMassless, kM, back, params, MTSE)

print Nstart
print backExitMinus

start_time = time.time()
# run solver for this triangle
t=np.linspace(Nstart,Nend, 1000)
threePt = MTSE.alphaEvolve(t,k1,k2,k3, backExitMinus,params,tols,True) # all data from three point run goes into threePt array
ppp= threePt[:,4+2*nF+6*2*nF*2*nF+1:]        

zzz= threePt[:,:5] # this the 3pt of zeta

fig3 = plt.figure(3)
# plots the reduced bispectrum
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')

print  time.time() - start_time
############################################################################################################################################




plt.show(fig1)
plt.show(fig2)
plt.show(fig3)


