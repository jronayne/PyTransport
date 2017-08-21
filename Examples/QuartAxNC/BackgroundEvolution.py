#################### generate power spectrum using PyTransQuartAxNC installed using the setup file which accompanies this one ####################
from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below
import timeit


location = "/home/jwr/Code//PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup
tols = np.array([10**-10,10**-10])

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransQuartAxNC as PyT  # import module as PyT (PyTransDQuad is quite long to type each time and it saves time to use a shorter name
                             # using a generic name PyT means the file can be more easily reused for a different example (once field values 
                             # etc are altered)
import PyTransScripts as PyS  # import the scripts module as PyS for convenience
###########################################################################################################################################
BB=10
fields = np.array([2.0,0.5-0.001])
nP=PyT.nP()   
nF=PyT.nF()   
params = np.zeros(nP)
nS = np.zeros(BB)
FNL=np.zeros([BB,1000])
zzOut=np.empty([BB,100])
zzzOut=np.empty([BB,100])
timess=np.empty([BB,100])
KO=np.empty([BB,100])
def get_last_non_zero_index(d, default=None):
		rev = (len(d) - idx for idx, item in enumerate(reversed(d), 1) if item)
		return next(rev, default)
#Sample over different values of the 2-sphere radius
for bb in range(BB):

	params[0]=1.*pow(10.,-10); params[1]=1.; params[2]=25.0**2.0*params[0]/4.0/math.pi**2;
	params[3]=9.3-0.01*bb

	V = PyT.V(fields,params) # calculate potential from some initial conditions
	dV=PyT.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)
	initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions using slow roll expression

	################################## run the background fiducial run #########################################################################
	Nstart = 0.0
	Nend = 65.0
	t=np.linspace(Nstart, Nend, 10000) # array at which output is returned
	back = PyT.backEvolve(t, initial, params, tols,True ) # The output is read into the back numpy array 

	numstart = (np.abs(back[:,0]-back[get_last_non_zero_index(back[:,0]),0]+64.0)).argmin()
	fieldsi=np.array([ back[numstart,1 ],back[numstart,2]])
	initiali = np.concatenate((fieldsi,np.array([0.,0.])))
	print numstart
	ti=np.linspace(Nstart, 65, 1000)
	back2 = PyT.backEvolve(ti, initiali, params, tols, False )
	
	fig1 = plt.figure(1)
	plt.plot(back2[:,0], back2[:,1 ], 'r',label=r'$\phi$',linewidth=2)
	plt.plot(back2[:,0],back2[:,2 ], 'g',label=r'$\chi$',linewidth=2)

plt.legend()
plt.title('Background Field Evolution')
plt.xlabel('N')
plt.ylabel('Field Value')
plt.show()