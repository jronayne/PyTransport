#################### Simple example of using PyTransCurve installed using the setup file which accompanies this one #########################

from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below
############################################################################################################################################

#This file contains simple examples of using the PyTransCurve
#It assumes the CurveSetup file has been run to install a double quadratic version of PyTransport
#It is recommended you restart the kernel before running this file to insure any updates to PyTransCurve are imported

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder folder
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this adds the other paths that PyTransport uses to the python path

import PyTransCurve as PyT  # import module as PyT 
                             # using a generic name PyT means the file can be more easily reused for a different example (once field values 
                             # etc are altered)
import PyTransScripts as PyS  # import the scripts module as PyS for convenience

###########################################################################################################################################
# Example 
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
tols = np.array([10**-10,10**-10])


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 28.0
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = PyT.backEvolve(t, initial, pvalue,tols, False) # The output is read into the back numpy array 

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

#k = PyS.kexitPhi(PhiExit, 1, back, pvalue, PyT) 
k = PyS.kexitN(10.0, back, pvalue, PyT) 

# other scales can then be defined wrt to k
##############################################################MB##############################################################################


################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds of massless sub-horizon evolution using a spline routine
NBMassless = 6.0
Nstart, backExitMinus = PyS.ICsBM(NBMassless, k, back, pvalue, PyT)  
print ("NStart= " + str(Nstart) )
print (k)
t=np.linspace(Nstart,Nend, 1000)  # array at which output is returned -- itial valye should correspond to initial field values
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an crude estimate for ns

twoPt = PyT.sigEvolve(t, k, backExitMinus,pvalue,tols,True) # puts information about the two point fuction in twoPt array
zz1=twoPt[:,1] # the second column is the 2pt of zeta
zz1a=zz1[-1] 
twoPt=PyT.sigEvolve(t, k+.1*k, backExitMinus,pvalue,tols,True) 
zz2=twoPt[:,1]
zz2a=zz2[-1]
n_s = (np.log(zz2a)-np.log(zz1a))/(np.log(k+.1*k)-np.log(k))+4.0
print (n_s)
fig2 = plt.figure(2)
plt.plot(t, zz1,'r')
yscale('log')
############################################################################################################################################

###################################### example bispectrum run ##############################################################################
# set three scales in FLS manner (using alpha beta notation)
alpha=0.
beta =1/3.

k1 = k/2 - beta*k/2.
k2 = k/4*(1+alpha+beta)
k3 = k/4*(1-alpha+beta)
#k1 =k; k2 =k ;k3 =k;
# find initial conditions for 4 e-folds of massless evolution for the largest k (which stays inside the longest)

kM = np.min(np.array([k1,k2,k3]))
Nstart, backExitMinus = PyS.ICsBM(NBMassless, kM, back, pvalue, PyT)

print (Nstart)
print (backExitMinus)

start_time = time.time()
# run solver for this triangle
t=np.linspace(Nstart,Nend, 1000)
t3p=np.linspace(Nstart,Nstart+.1, 2)
threePt = PyT.alphaEvolve(t,k1,k2,k3, backExitMinus,pvalue,tols,True) # all data from three point run goes into threePt array
alpha= threePt[:,5+2*nF+6*2*nF*2*nF:]        
zzz= threePt[:,:5] # this the 3pt of zeta

fig4 = plt.figure(4)
for ii in range(0,2):
    for jj in range(0,2):
        for kk in range(0,2):
            plt.plot(t, np.abs(alpha[:,ii + 2*nF*jj + 2*nF*2*nF*kk]))
title(r'$\alpha$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Absolute 3pt field correlations', fontsize=20);
xlabel(r'$N$', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15);

fig3 = plt.figure(3)
# plots the reduced bispectrum
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')

print  (time.time() - start_time)
############################################################################################################################################


plt.show(fig1)
plt.show(fig2)
plt.show(fig3)


