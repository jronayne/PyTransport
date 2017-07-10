####################################### Setup file for the Non-Canonical quasi-single field Model ######################################################

import sympy as sym
import numpy as np
import math
import sys
from pylab import *

from gravipy import *

############################################################################################################################################

location = "/home/jwr/Code/PyTransport/" # this should be the location of the PyTransport folder
sys.path.append(location)  # we add this location to the python path

import PyTransSetup as PySet

### Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ###
### Restart the python kernel after running this file
nF=2
nP=4

f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
G=Matrix( [[ 1.0,0], [0,f[0]**2.0] ] )

V = 10.0**(-10.0) * (1.0 + 29.0/120. *math.pi * f[1] +  1./2. * p[1] * (f[0] - p[0])**2 +  1./3./2. * p[2] * (f[0] - p[0])**3 +  1./4./3./2. * p[3] * (f[0] - p[0])**4)
                  
PySet.tol(1e-8,1e-8)
#The last argument is for whether the sympy's simplify is used to on derivatives of the potential and field geometric quantites.
#Caution is recomended as sympy's simplify is known to have bugs. Simplification can increase the speed of numerically evolutions, but at the cost of compling more slowly.
PySet.potential(V,nF,nP,False,G) # writes this potential into c file when run
##Set second argument to True to use the non-canonical set-up
PySet.compileName("CurveNC",True) # this compiles the module with the new potential and places it in the location folder, and adds this folder to the path ready for use
############################################################################################################################################

