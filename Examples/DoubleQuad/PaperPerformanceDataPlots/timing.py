#################### generate power spectrum using MTeasyDQuad installed using the setup file which accompanies this one #########################

from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import timeit  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below
############################################################################################################################################

#This file contains simple examples of using the MTeasyDQuad 
#It assumes the DQuadSetup file has been run to install a double quadratic version of MTeasyPy
#It is recommended you restart the kernel before running this file to insure any updates to MTeasyPyDQuad are imported

location = "/Users/david/Dropbox/MTeasyDist/MTeasy/" # this should be the location of the MTeasy folder 
sys.path.append(location) # sets up python path to give access to MTeasySetup

import MTeasySetup
MTeasySetup.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyDQuad as MTE  # import module as MTE (MTeasyDQuad is quite long to type each time and it saves time to use a shorter name
                             # using a generic name MTE means the file can be more easily reused for a different example (once field values 
                             # etc are altered)
import MTeasyScripts as MTS  # import the scripts module as MTS for convenience

###########################################################################################################################################

########################### set up initial conditions for background run ##################################################################
fields = np.array([12.9, 10.0]) # we set up a numpy array which contains the values of the fields

nP=MTE.nP()   # the .nP() function gets the number of parameters needed for the potential -- this can be used as a useful crosscheck
params = np.zeros(nP) 
params[0]=10.0**(-5.0); params[1]=9.0*10.0**(-5.0) # we set up numpy array which contains values of the parameters

nF=MTE.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)

V = MTE.V(fields,params) # calculate potential from some initial conditions
dV=MTE.dV(fields,params) # calculate derivatives of potential 

initial = np.concatenate((fields, -dV/np.sqrt(3* V))) # sets an array containing field values and there derivative in cosmic time 
                                                      # (set using the slow roll equation)
############################################################################################################################################



################################## run the background fiducial run ################################################################
Nstart = 0.0
Nend = 65.0
t=np.linspace(Nstart, Nend, 1000)

back = MTE.backEvolve(t, initial, params)
###########################################################################################################################################


###########################################################################################################################################
# set up k values 
kt = MTS.kexitN(24.5,back,params,MTE)


alpha = 0.0; beta = .99
kt=3.*kt
k1 = kt/2. - beta*kt/2.; k2 = kt/4.*(1.+alpha+beta); k3 = kt/4.*(1.-alpha+beta)
kM = min(k1,k2,k3)
NBefore = np.linspace(2,8,20)
zzzOut = np.zeros(np.size(NBefore)); zz1Out = np.zeros(np.size(NBefore)); zz2Out = np.zeros(np.size(NBefore)); zz3Out = np.zeros(np.size(NBefore)); times = np.zeros(np.size(NBefore));


for ii in range(0,np.size(NBefore)):
    print ii
    Nstart, backExitMinus = MTS.ICsBE(NBefore[ii], kM, back, params, MTE)
    timebefore = timeit.default_timer()
    t=np.linspace(Nstart,Nend, 10)  # array at which output is returned -- initial value should correspond to initial field values
    threePt = MTE.alphaEvolve(t,k1,k2,k3, backExitMinus,params,0)
    time = timeit.default_timer() - timebefore
    zzzOut[ii] = threePt[-1,4]; zz1Out[ii] = threePt[-1,1]; zz2Out[ii] = threePt[-1,2]; zz3Out[ii] = threePt[-1,3];
    times[ii] = time
  
fig1 = plt.figure(1)
plt.scatter(NBefore,times)
#ylim([9./10.*np.min(times),10./9.*np.max(times)])
grid(True); plt.legend(fontsize=15); ylabel(r'time taken', fontsize=20); 
xlabel(r'e-folds before exit', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15); plt.savefig("times4.png")
  
fig2 = plt.figure(2)
DBi=zzzOut*k1**2*k2**2*k3**2
plt.scatter(NBefore,DBi)
ylim([9./10.*np.min(DBi),10./9.*np.max(DBi)])
grid(True); plt.legend(fontsize=15); ylabel(r'dimensionless bisecptrum', fontsize=20); 
xlabel(r'e-folds before exit', fontsize=15); grid(True); plt.legend(fontsize=15); plt.savefig("DBiCon4.png")
    
  
np.savetxt('data/timesSq2.dat', (zzzOut,zz1Out,zz2Out,zz3Out,times))   # x,y,z equal sized 1D arrays

