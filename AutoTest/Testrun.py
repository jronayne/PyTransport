#Test file for use with nose.
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below

location = "/home/jwr/Code/nosetest/PyTransport-master/PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup
tols = np.array([10**-10,10**-10])

import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTransport uses
import PyTransQuartAxNC as PyT  
import PyTransScripts as PyS  # import the scripts module as PyS for convenience

fields = np.array([2.0,0.5-0.001])
nP=PyT.nP()   
nF=PyT.nF()   
params = np.zeros(nP)

params[0]=1.*pow(10.,-10); params[1]=1.; params[2]=25.0**2.0*params[0]/4.0/math.pi**2;params[3]=9.1
V = PyT.V(fields,params) # calculate potential from some initial conditions
dV=PyT.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions using slow roll expression
Nstart = 0.0
Nend = 65.0
t=np.linspace(Nstart, Nend, 10000) # array at which output is returned
back = PyT.backEvolve(t, initial, params, tols,True) # The output is read into the back numpy array 

##test 1
def test_Background():
	BR=back[:,1:]
	BD=np.loadtxt('TEST_backevo.out')
	Run1=np.array([str(BR[i,k])[0:10] for i in range(9862) for k in range(4)])
	Data1=np.array([str(BD[i,k])[0:10] for i in range(9862) for k in range(4)])
	assert np.array_equal(Run1,Data1)


k = PyS.kexitN(14.0, back, params, PyT)
NB = 5.0
Nstart, backExitMinus = PyS.ICsBE(NB, k, back, params, PyT) #find conditions for 5 e-folds before horizon crossing of k mode
tsig=np.linspace(Nstart,Nend, 1000)
twoPt = PyT.sigEvolve(tsig, k, backExitMinus,params,tols, True) # puts information about the two point fuction in twoPt array
zz1=twoPt[:,1] # the second column is the 2pt of zeta
sigma = twoPt[:,1+1+2*nF:] # the last 2nF* 2nF columns correspond to the evolution of the sigma matrix
zz1a=zz1[-1] # the value fo the power spectrum for this k value at the end of the run
twoPt=PyT.sigEvolve(tsig, k+.1*k, backExitMinus,params,tols, True)
zz2=twoPt[:,1]
zz2a=zz2[-1]
n_s = (np.log(zz2a)-np.log(zz1a))/(np.log(k+.1*k)-np.log(k))+4.0

##test 2
def test_2pt():
	ZZR=(zz1,zz2)
	ZZD=np.loadtxt('TEST_zz.out')
	Run2=np.array([str(ZZR[k][i])[0:10] for i in range(1000) for k in range(2)])
	Data2=np.array([str(ZZD[k][i])[0:10] for i in range(1000) for k in range(2)])
	assert np.array_equal(Run2,Data2)


alpha=0.0
beta =1/3.
k1 = k/2 - beta*k/2. ; k2 = k/4*(1+alpha+beta) ; k3 = k/4*(1-alpha+beta)
kM = np.min(np.array([k1,k2,k3]))
Nstart, backExitMinus = PyS.ICsBE(NB, kM, back, params, PyT)

talp=np.linspace(Nstart,Nend, 1000)
threePt = PyT.alphaEvolve(talp,k1,k2,k3, backExitMinus,params,tols, True)
alpha= threePt[:,1+4+2*nF+6*2*nF*2*nF:]        # this now contains the 3pt of the fields and field derivative pertruabtions
zzz= threePt[:,1:5]
fnl = 5.0/6.0*zzz[:,3]/(zzz[:,1]*zzz[:,2]  + zzz[:,0]*zzz[:,1] + zzz[:,0]*zzz[:,2])

#test 3
def test_3pt():
	ZZZR=(zzz)
	ZZZD=np.loadtxt('TEST_zzz.out')
	Run3=np.array([str(ZZZR[i,k])[0:10] for i in range(1000) for k in range(4)])
	Data3=np.array([str(ZZZD[i,k])[0:10] for i in range(1000) for k in range(4)])
	assert np.array_equal(Run3,Data3)