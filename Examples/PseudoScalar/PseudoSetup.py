####################################### Setup file for the PyTransPseudo Model ######################################################

import sympy as sym    # we import the sympy package
import math            # we import the math package (not used here, but has useful constants such math.pi which might be needed in other cases)    
import sys             # import the sys module used below
from pylab import *
from gravipy import *

############################################################################################################################################

# if using an integrated environment we recommend restarting the python console after running this script to make sure updates are found 

location = "/home/jwr/Code/PyTransport/" # this should be the location of the PyTransport folder 
sys.path.append(location)  # we add this location to the python path

import PyTransSetup  # the above commands allows python to find the PyTransSetup module and import it

############################################################################################################################################

nF=3  # number of fields needed to define the PseudoScalar potential
nP=3  # number of parameters needed to define the PseudoScalar potential
f=sym.symarray('f',nF)   # an array representing the nF fields present for this model
p=sym.symarray('p',nP)   # an array representing the nP parameters needed to define this model (that we might wish to change) if we don't 
V= p[0] * f[0]**2 + p[1] * f[1]**2 + p[2] * f[2]**2 # this is the potential written in sympy notation
R=(0.9)/((sym.cosh(2.0*((1.0*f[0])-7.0)/0.12))**(2.0))
G=Matrix([[1,R,0],[R,1,0],[0,0,1]]) # selecting the field space metric in this instance.
PyTransSetup.tol(1e-12,1e-12)   # set tols for the numerical integration of the 3 point function (if not run this will remain as it the
                               #last time set)

PyTransSetup.potential(V,nF,nP,False,G) # writes this potential and its derivatives into C++ file potential.h when run

PyTransSetup.compileName("Pseudo",True) # this compiles a python module using the C++ code, including the edited potential.h file, called PyTransPseudo
                                 # and places it in the location folder, ready for use

############################################################################################################################################