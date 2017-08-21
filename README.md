# PyTransport
PyTransport release 2.0 (2017).
This code has been written by Dr. David J. Mulryne and John W. Ronayne.

This upload contains the PyTransport code, as well as examples and many of the files and figures generated during testing. 

PyTransport constitutes a straightforward code written in C++  together with Python scripts which automatically edit, compile and run the C++ code as a Python module. The code is intended to be a reusable resource for inflationary cosmology. It enables users to quickly create a complied Python module(s) for any given model(s) of multi-field inflation. The primary function of the complied module is to calculate the power-spectrum and bi-spectrum of inflationary perturbations produced by multi-field inflation. 

# Update in Version 2.0 
New modules have been added and modifications to existing functions have been made to support non-canonical models of Inflation. In addition, new arguments in some functions are available to the end user which allow more control of the speed and output of the simulations. These new features modify the input of existing evolution functions to allow; the ability of the integrator to run until the end of inflation, for tolerances to be set for each individual function and the simplification of the model potential and field-metric quantities to speed up numerical evaluation at the sacrifice of a slower compilation.
In addition error feedback has been improved and storage of model quantities has been made more efficient. 

PyTransport has been developed on OS X and Linux Ubuntu using Python 2.7, and is intended for use on Unix based systems.

# Licencing #
PyTransport is distributed under the GNU General Public License version 3, or (at your option) any later version. This license is bundled with the source code as LICENSE.txt.
