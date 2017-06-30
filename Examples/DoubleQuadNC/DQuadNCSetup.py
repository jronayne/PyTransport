####################################### Setup file for the double quadratic potential with 2-sphere field space metric######################################################

import sympy as sym    # we import the sympy package
import math            # we import the math package (not used here, but has useful constants such math.pi which might be needed in other cases)    
import sys             # import the sys module used below
from gravipy import *  # we import the gravipy package to write the field metric
############################################################################################################################################

# if using an integrated environment we recommend restarting the python console after running this script to make sure updates are found 

location = "/home/jwr/Code/PYT4/PyTransport/" # this should be the location of the PyTrans folder
sys.path.append(location)  # we add this location to the python path

import PyTransSetup  # the above commands allows python to find the PyTransSetup module and import it

############################################################################################################################################

nF=2  # number of fields needed to define the double quadratic potential
nP=3  # number of parameters needed to define the double quadtartic potential
f=sym.symarray('f',nF)   # an array representing the nF fields present for this model
p=sym.symarray('p',nP)   # an array representing the nP parameters needed to define this model (that we might wish to change) if we don't 
                         # wish to change them they could be typed explicitly in the potential below

V = 1./2. * p[0]**2.0 *f [0]**2.0  +  1./2. * p[1]**2.0 * f[1]**2.0   # this is the potential written in sympy notation
G=Matrix( [[p[2]**2.0, 0], [0, p[2]**2.0*sym.sin(f[0])**2.0] ] ) # this is the field metric written in sympy notation
#PyTransSetup.tol(1e-8,1e-8)   # set tols for the numerical integration of the 3 point function (if not run this will remain as it the
                               #last time set)
#The last argument is for whether the sympy's simplify is used to on derivatives of the potential and field geometric quantites.
#Caution is recomended as sympy's simplify is known to have bugs. Simplification can increase the speed of numerically evolutions, but at the cost of compling more slowly.
PyTransSetup.potential(V,nF,nP,G) # writes this potential and its derivatives into C++ file potential.h when run

PyTransSetup.compileName("DQuadNC",True) # this compiles a python module using the C++ code, including the edited potential.h file, called PyTransDQuad
                                 # and places it in the location folder, ready for use

############################################################################################################################################

