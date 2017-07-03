#################### generate power spectrum using PyTransCurveNC installed using the setup file which accompanies this one ####################
from matplotlib import pyplot as plt   # import package for plotting
from pylab import *  # contains some useful stuff for plotting
import time  # imports a package that allows us to see how long processes take
import math  # imports math package
import numpy as np # imports numpu package as np for short
import sys  # imports sys package for sue below


############################################################################################################################################


############################################################################################################################################
#This file contains simple examples of using the PyTransCurveNC 
#It assumes the PyTransCurveNC file has been run to install a non-canonical curve version of PyTrans
#It is recommended you restart the kernel before running this file to insure any updates to PyTransCurveNC are imported

location = "/home/jwr/Code/PYT4/PyTransport/" # this should be the location of the PyTrans folder 
sys.path.append(location) # sets up python path to give access to PyTransSetup
tols = np.array([10**-12,10**-12])
import PyTransSetup
PyTransSetup.pathSet()  # this sets the other paths that PyTrans uses

import PyTransCurveNC as PyT  # import module as PyT (PyTransCurveNC is quite long to type each time and it saves time to use a shorter name
                             # using a generic name PyT means the file can be more easily reused for a different example (once field values 
                             # etc are altered)
import PyTransScripts as PyS  # import the scripts module as PyS for convenience
###########################################################################################################################################

omega = pi/30.0
R0 = np.sqrt(10.0**(-10)/3.0) / (omega *np.sqrt(10.0**(-9)))
fields = np.array([np.sqrt(R0**2.0+((1e-2)*R0)**2.0),math.atan2((1e-2)*R0,-R0)])#  (1e-0)*R0]) # we set up a numpy array which contains the values of the fields
nP=PyT.nP()   # the .np function gets the number of parameters needed for the potential -- this can be used as a useful crosscheck
pvalue = np.zeros(nP)
pvalue[1]=1./np.sqrt(3.0); pvalue[2]=1.0/10.0**(-5); pvalue[3]=1.0/omega**.5 * 1./2. * (10.0**(-10))**(-3./4.)
pvalue[0]=R0

nF=PyT.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)

V = PyT.V(fields,pvalue) # calculate potential from some initial conditions
dV = PyT.dV(fields,pvalue) # calculate derivatives of potential

initial = np.concatenate((fields, np.array([0,0]))) # sets an array containing field values and there derivative in cosmic time
# (set using the slow roll equation)

############################################################################################################################################
################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 31
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = PyT.backEvolve(t, initial, pvalue,tols, True) # The output is read into the back numpy array
np.savetxt('back.out',back)
#plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r',label=r'$\theta$')
plt.plot(back[:,0], back[:,1 ], 'g',label=r'$R$')
title(r'Background field evolution in polar frame',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Field value', fontsize=15); 
xlabel(r'$N$', fontsize=15); grid(True);  plt.legend(fontsize=15);
fig2 = figure(2)
plt.plot(back[:,0], back[:,1 ]*np.cos(back[:,2 ]), 'r',label=r'$Y$')
plt.plot(back[:,0], back[:,1 ]*np.sin(back[:,2 ]), 'g',label=r'$X$')
title(r'Background field evolution in cartesian frame',fontsize=15); grid(True); plt.legend(fontsize=15,loc='upper left'); ylabel(r'Field value', fontsize=15); 
xlabel(r'$N$', fontsize=15); grid(True);  plt.legend(fontsize=15,loc='upper left');
fig3 = figure(3)
plt.plot(back[:,1 ]*np.sin(back[:,2 ]), back[:,1 ]*np.cos(back[:,2 ]), 'or')
############################################################################################################################################
k = PyS.kexitN(10.0, back, pvalue, PyT) 
plt.show()
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

tols = np.array([10**-8,10**-8])
twoPt = PyT.sigEvolve(t, k, backExitMinus,pvalue,tols, True) # puts information about the two point fuction in twoPt array
zz1=twoPt[:,1] # the second column is the 2pt of zeta
zz1a=zz1[-1] 
twoPt=PyT.sigEvolve(t, k+.1*k, backExitMinus,pvalue,tols, True) 
zz2=twoPt[:,1]
zz2a=zz2[-1]
n_s = (np.log(zz2a)-np.log(zz1a))/(np.log(k+.1*k)-np.log(k))+4.0
np.savetxt('sigma.out',twoPt)
print (n_s)
fig4 = plt.figure(4)
plt.plot(t, zz1,'r')
title(r'$<\zeta\zeta>$ Two-point function',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Correlation amplitude', fontsize=15); 
xlabel(r'$N$', fontsize=15); grid(True);  plt.legend(fontsize=15);
yscale('log')

sigma = twoPt[:,1+1+2*nF:]
fig10=plt.figure(10)
for ii in range(0,2):
    for jj in range(0,2):
        plt.plot(twoPt[:,0], np.abs(sigma[:,ii + 2*nF*jj]))
title(r'$\Sigma$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Aboslute 2pt field correlations', fontsize=20); 
xlabel(r'$N$', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15);
np.savetxt('sigma.out',sigma)
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

threePt = PyT.alphaEvolve(t,k1,k2,k3, backExitMinus,pvalue,tols, True) # all data from three point run goes into threePt array

alpha= threePt[:,5+2*nF+6*2*nF*2*nF:]        
zzz= threePt[:,:5] # this the 3pt of zeta
zzz1= threePt[:,4]
fig4 = plt.figure(9)
plt.plot(np.abs(zzz1))
yscale('log')
plt.show()
np.savetxt('zzz.out',zzz1)
np.savetxt('alpha.out',alpha)

fig4 = plt.figure(4)
for ii in range(0,2):
    for jj in range(0,2):
        for kk in range(0,2):
            plt.plot(t, np.abs(alpha[:,ii + 2*nF*jj + 2*nF*2*nF*kk]))
title(r'$\alpha$ evolution',fontsize=15); grid(True); plt.legend(fontsize=15); ylabel(r'Absolute 3pt field correlations', fontsize=20);
xlabel(r'$N$', fontsize=15); grid(True); yscale('log'); plt.legend(fontsize=15);

	
fig5 = plt.figure(5)
# plots the reduced bispectrum
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')

print  (time.time() - start_time)
############################################################################################################################################


plt.show()

############################################################################################################################################
