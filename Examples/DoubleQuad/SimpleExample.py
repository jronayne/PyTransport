#################### Simple example of using PyTransDQuad installed using the setup file which accompanies this one #########################

from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below
############################################################################################################################################

#This file contains simple examples of using the PyTransDQuad 
#It assumes the DQuadSetup file has been run to install a double quadratic version of PyTrans
#It is recommended you restart the kernel before running this file to insure any updates to PyTransDQuad are imported

location = "/home/jwr/Code/PyTransport/" # this should be the location of the PyTransport folder folder
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # this adds the other paths that PyTrans uses to the python path

import PyTransDQuad as PyT  # import module as PyT (PyTransDQuad is quite long to type each time and it saves time to use a shorter name
                             # using a generic name PyT means the file can be more easily reused for a different example (once field values 
                             # etc are altered)
import PyTransScripts as PyS  # import the scripts module as PyS for convenience

###########################################################################################################################################



# Example 

########################### set some field values and field derivatives in cosmic time ####################################################

fields = np.array([12.0, 12.0]) # we set up a numpy array which contains the values of the fields

nP=PyT.nP()   # the .np function gets the number of parameters needed for the potential -- this can be used as a useful crosscheck
pvalue = np.zeros(nP) 
pvalue[0]=10.0**(-5.0); pvalue[1]=9.0*10.0**(-5) # we set up numpy array which contains values of the parameters

nF=PyT.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)

V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV=PyT.dV(fields,pvalue) # calculate derivatives of potential

initial = np.concatenate((fields, -dV/np.sqrt(3* V))) # sets an array containing field values and there derivative in cosmic time 
                                                      # (set using the slow roll equation)

############################################################################################################################################



################################## run and plot the background fiducial run ################################################################
Nstart = 0.0
Nend = 80.0
t=np.linspace(Nstart, Nend, 1000)

tols = np.array([10**-8,10**-8])
back = PyT.backEvolve(t, initial, pvalue,tols, True)


fig1=plt.figure(1)
plt.plot(back[:,0], back[:,1], 'g')
plt.plot(back[:,0], back[:,2], 'r')
title(r'Background evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Fields', fontsize=20); xlabel(r'$N$', fontsize=15)
grid(True); plt.legend(fontsize=15); plt.savefig("DQ1.png")
############################################################################################################################################


print(back[-1,1])
print(size(back[:,1]))
############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t

NExit = 15.0
k = PyS.kexitN(NExit, back, pvalue, PyT) 

# other scales can then be defined wrt to k
############################################################################################################################################


################################# example 2pt run ##########################################################################################

NB = 6.0
Nstart, backExitMinus = PyS.ICsBE(NB, k, back, pvalue, PyT) #find conditions for 6 e-folds before horizon crossing of k mode


tsig=np.linspace(Nstart,back[-1,0], 1000)  # array of times (e-folds) at which output is returned -- initial time should correspond to initial field
                                     # and velocity values which will be fed in to the functions which evolve correlations


# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be
# repeated to build up the spectrum, here we run twice to get an crude estimate for ns
twoPt = PyT.sigEvolve(tsig, k, backExitMinus,pvalue,tols, True) # puts information about the two point fuction in twoPt array
zz1=twoPt[:,1] # the second column is the 2pt of zeta
sigma = twoPt[:,1+1+2*nF:] # the last 2nF* 2nF columns correspond to the evolution of the sigma matrix
zz1a=zz1[-1] # the value of the power spectrum for this k value at the end of the run

twoPt=PyT.sigEvolve(tsig, k+.1*k, backExitMinus, pvalue,tols, True)
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
beta =1.0/3.0

k1 = k/2 - beta*k/2. ; k2 = k/4*(1+alpha+beta) ; k3 = k/4*(1-alpha+beta)


# find initial conditions for 6 e-folds before the smallest k (which exits the horizon first) crosses the horizon
kM = np.min(np.array([k1,k2,k3]))
Nstart, backExitMinus = PyS.ICsBM(NB, kM, back, pvalue, PyT)


# run the three point evolution for this triangle
talp=np.linspace(Nstart,back[-1,0], 1000)
threePt = PyT.alphaEvolve(talp,k1,k2,k3, backExitMinus,pvalue,tols,True) # all data from three point run goes into threePt array
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
