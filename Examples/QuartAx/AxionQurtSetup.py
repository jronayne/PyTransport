####################################### Setup file for the axion quartic example ##########################################################
import sympy as sym
import numpy as np
import math
import sys
############################################################################################################################################

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTransport folder
sys.path.append(location)  # we add this location to the python path

import PyTransSetup

### Sets potential and compiles PyTransport, users may prefer to do this only once in a separate file (or comment after running below once) ###
### Restart the python kernel after running this file
nF=2
nP=3
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./4. * p[0] * f[0]**4 + p[2] * (1-sym.cos(2*math.pi * f[1] / p[1]))


PyTransSetup.tol(1e-8,1e-8)
PyTransSetup.potential(V,nF,nP) # writes this potential into c file when run

PyTransSetup.compileName("AxQrt") # this compiles the module with the new potential and places it in the location folder, and adds this folder to the path ready for use
############################################################################################################################################

