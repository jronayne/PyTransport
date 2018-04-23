#################### Simple example of using PyTransAxQrt installed using the setup file which accompanies this one #########################

from matplotlib import pyplot as plt
from pylab import *  # contains some useful stuff for plotting
import timeit
import math 
import numpy as np
import sys 
############################################################################################################################################

#This file contains simple examples of using the PyTransport package for the heavy field example of Langlois.
#It assumes the PyTransAxQrt file has been run to install a AxQrt version of PyTransPy
#It is recommended you restart the kernel to insure any updates to PyTransAxQrt are imported 

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # his add sets the other paths that PyTransport uses

import PyTransAxQrt as PyT;  # import module
import PyTransScripts as PyS;

# Example

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
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
############################################################################################################################################


############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
#PhiExit = -100.*math.sqrt(6.) + shift + .85

#k = PyS.kexitPhi(PhiExit, 1, back, params, MTLH) 
k = PyS.kexitN(14.0, back, params, PyT)

# other scales can then be defined wrt to k
############################################################################################################################################


################################# example 2pt run ##########################################################################################

NB = 4.5
Nstart, backExitMinus = PyS.ICsBE(NB, k, back, params, PyT) #find conditions for 5 e-folds before horizon crossing of k mode


tsig=np.linspace(Nstart,Nend, 1000)  # array at which output is returned -- initial value should correspond to initial field values


# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an crude estimate for ns
twoPt = PyT.sigEvolve(tsig, k, backExitMinus,params, tols,True) # puts information about the two point fuction in twoPt array
zz1=twoPt[:,1] # the second column is the 2pt of zeta
sigma = twoPt[:,1+1+2*nF:] # the last 2nF* 2nF columns correspond to the evolution of the sigma matrix
zz1a=zz1[-1] # the value fo the power spectrum for this k value at the end of the run

twoPt=PyT.sigEvolve(tsig, k+.1*k, backExitMinus,params, tols,True)
zz2=twoPt[:,1]
zz2a=zz2[-1]
n_s = (np.log(zz2a)-np.log(zz1a))/(np.log(k+.1*k)-np.log(k))+4.0


fig2=plt.figure(2)
for ii in range(0,2):
    for jj in range(0,2):
        plt.plot(twoPt[:,0], np.abs(sigma[:,ii + 2*nF*jj]))
title(r'$\Sigma$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Aboslute 2pt field correlations', fontsize=20); 
xlabel(r'$N$', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15); plt.savefig("DQ2.png")


fig3=plt.figure(3)
plt.plot(tsig, zz1[:])
title(r'$P_\zeta$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'$P_\zeta(k)$', fontsize=20); 
xlabel(r'$N$', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15); plt.savefig("DQ3.png")

#############################################################################################################################################


###################################### example bispectrum run ##############################################################################

# set three scales in FLS manner (using alpha, beta notation)
alpha=0.0
beta =1/3.

k1 = k/2 - beta*k/2. ; k2 = k/4*(1+alpha+beta) ; k3 = k/4*(1-alpha+beta)


# find initial conditions for 5 e-folds before the largest k (which stays inside the longest) crosses the horizon
kM = np.min(np.array([k1,k2,k3]))
Nstart, backExitMinus = PyS.ICsBE(NB, kM, back, params, PyT)


# run the three point evolution for this triangle
talp=np.linspace(Nstart,Nend, 1000)
timebefore = timeit.default_timer()
threePt = PyT.alphaEvolve(talp,k1,k2,k3, backExitMinus,params,1) # all data from three point run goes into threePt array
time = timeit.default_timer()-timebefore
print time
alpha= threePt[:,1+4+2*nF+6*2*nF*2*nF:]        # this now contains the 3pt of the fields and field derivative pertruabtions
zzz= threePt[:,1:5] # this contains the evolution of two point of zeta for each k mode involved and the 3pt of zeta

fig4 = plt.figure(4)
for ii in range(0,2):
    for jj in range(0,2):
        for kk in range(0,2):
            plt.plot(talp, np.abs(alpha[:,ii + 2*nF*jj + 2*nF*2*nF*kk]))
title(r'$\alpha$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Absolute 3pt field correlations', fontsize=20); 
xlabel(r'$N$', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15); plt.savefig("DQ4.png")

fig5=plt.figure(5)
fnl = 5.0/6.0*zzz[:,3]/(zzz[:,1]*zzz[:,2]  + zzz[:,0]*zzz[:,1] + zzz[:,0]*zzz[:,2])
plt.plot(talp, fnl,'r')
title(r'$f_{NL}$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'$f_{NL}$', fontsize=20); 
xlabel(r'$N$', fontsize=15); grid(True); plt.legend(fontsize=15); plt.savefig("DQ5.png")


############################################################################################################################################




plt.show(fig1);plt.show(fig2);plt.show(fig3);plt.show(fig4);plt.show(fig5)
