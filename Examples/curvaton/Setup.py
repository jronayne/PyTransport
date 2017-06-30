####################################### Setup file for the Step Potential example of Chen et al. ###########################################
import sympy as sym
import numpy as np
import math
import sys
############################################################################################################################################

location = "/Users/david/Dropbox/PyTransportDist/PyTransport/" # this should be the location of the PyTrans folder
sys.path.append(location)  # we add this location to the python path

import PyTransSetup

### Sets potential and compiles PyTrans, users may prefer to do this only once in a separate file (or comment after running below once) ###
nF=2
nP=3
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)

## example step
V = 1.0/2.0 *p[0]**2*f[0]**2 + 1.0/2.0 *p[1]**2*f[1]**2 + 1.0/4.0 *p[2]*f[1]**4

PyTransSetup.tol(1e-8,1e-8)

PyTransSetup.potential(V,nF,nP) # differentiates this potential and writes this potential and derivatives into c file when run (can be a 
                               # little slow, and so one may not wish to run if recompiling to alater other properties such as tols) 

PyTransSetup.compileName("Cton") # this compiles the module with the new potential and places it in the location folder, and adds this folder to the path ready for use
############################################################################################################################################

