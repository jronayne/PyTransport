####################################### Setup file for the quasi-single field example  ####################################################
import sympy as sym
import numpy as np
import math
import sys
############################################################################################################################################

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder
sys.path.append(location)  # we add this location to the python path

import PyTransSetup as PySet

### Sets potential and compiles PyTransport, users may prefer to do this only once in a separate file (or comment after running below once) ###
### Restart the python kernel after running this file

nF=2
nP=4
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)

#V = 10.0**(-10.0) * (10.0 - (2.0 * p[0])**(1./2.) * sym.atan(f[0]/f[1]) +  1./2. * p[1]**2.0*(4.0-4.0*(f[1]**2.0 + f[0]**2.0)**(1/2.) + (f[1]**2.0 + f[0]**2.0)))
V = 10.0**(-10.0) * (1.0 + 29.0/120. *math.pi * sym.atan2(f[1],f[0])
                     +  1./2. * p[1] * ((f[1]**2 + f[0]**2)**(1./2.) - p[0])**2
                     +  1./3./2. * p[2] * ((f[1]**2 + f[0]**2)**(1./2.) - p[0])**3
                     +  1./4./3./2. * p[3] * ((f[1]**2 + f[0]**2)**(1./2.) - p[0])**4)
                     
PySet.potential(V,nF,nP) # writes this potential into c file when run

PySet.compileName("Curve") # this compiles the module with the new potential and places it in the location folder, and adds this folder to the path ready for use
############################################################################################################################################

