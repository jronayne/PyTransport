from matplotlib import pyplot as plt
from pylab import *  # contains some useful stuff for plotting
import timeit
import math 
import numpy as np
import sys 
############################################################################################################################################

#This file contains simple examples of using the PyTrans package for the heavy field example of Langlois.
#It assumes the LangHeavySeptup file has been run to install a LH version of PyTransPy
#It is recommended you restart the kernel to insure any updates to PyTransPyLH are imported 

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder
sys.path.append(location) # sets up python path to give access to PyTransSetup

import PyTransSetup
PyTransSetup.pathSet()  # his add sets the other paths that PyTrans uses

import PyTransCton as PyT;  # import module
import PyTransScripts as PyS;

# Example

########################### set initial field values and parameters for a simple example run ###################################################
nF=PyT.nF() # gets number of fields (useful check)
nP=PyT.nP() # gets number of parameters needed (useful check)

fields = np.array([5.0,0.01])
tols = np.array([10**-10,10**-10])

params = np.zeros(nP)
params[0]=1.; params[1]=.01*10.0**(-2); params[2]=(.5*10.0**(-2))**2 / fields[1]**3;

V = PyT.V(fields,params) # calculate potential from some initial conditions
dV=PyT.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions using slow roll expression
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 17
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = PyT.backEvolve(t, initial, params,tols, False) # The output is read into the back numpy array 
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
############################################################################################################################################
plt.show()