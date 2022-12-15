# PyTransport
PyTransport release 2.0 (2017).
This code has been written by Dr. David J. Mulryne and John W. Ronayne.

This upload contains the PyTransport code, as well as examples and many of the files and figures generated during testing.  

PyTransport constitutes a straightforward code written in C++  together with Python scripts which automatically edit, compile and run the C++ code as a Python module. The code is intended to be a reusable resource for inflationary cosmology. It enables users to quickly create a complied Python module(s) for any given model(s) of multi-field inflation. The primary function of the complied module is to calculate the power-spectrum and bi-spectrum of inflationary perturbations produced by multi-field inflation. 

# Update in Version 2.0 
New modules have been added and modifications to existing functions have been made to support non-canonical models of Inflation. In addition, new arguments in some functions are available to the end user which allow more control of the speed and output of the simulations. These new features modify the input of existing evolution functions to allow; the ability of the integrator to run until the end of inflation, for tolerances to be set for each individual function and the simplification of the model potential and field-metric quantities to speed up numerical evaluation at the sacrifice of a slower compilation.
In addition error feedback has been improved and storage of model quantities has been made more efficient. 

PyTransport has been developed on OS X and Linux Ubuntu using Python 2.7, and is intended for use on Unix based systems.
# Quick installation Guide
In order to compile PyTransport some prerequisites are required:
* A working Python installation (We recommend Python 2.7 however we have attempted to ensure compatibility with versions of Python 3).
* The following python packages, Numpy, Matplotlib, SciPy, Gravipy (v 0.1.0 *Note: compatiability issues with v0.2.0*), SymPy, Distutils, Math and Sys.
The simplist way to install these packages is by using pip, e.g.
```sh
pip install numpy
```
* Optional packages: Mpi4Py (also requiring [openMPI](https://wiki.helsinki.fi/display/HUGG/Open+MPI+install+on+Mac+OS+X) for distributed computing and Mayavi for 3D Bispectra plots.
* C++ compiler.

## Installing PyTransport
Download the repository and move it to a convenient location on your computer's file system.
Each model requires a separate installation. The `PyTransport-master/PyTransport/Examples/` folder contains subfolders of sample inflationary models, within them are the `ModelSetup.py` scripts.
Each setup scripts can be modified for a particular model (for a full description on modifying the setup file see the user guide in the `PyTransport-master/Docs/` folder).

For example, let's say you want to install a PyTransport module for the non-canonical Quartic-Axion model. 
First, from the shell, navigate into the folder `PyTransport-master/PyTransport/Examples/QuartAxNC/` and open the file `QuartAxNCsetup.py`.

You will first need to modify the following line which specifies the location of your PyTransport folder on your system.
```python
location = "/path/to/PyTransport/" # this should be the location of the PyTransport folder 
```
Note, if you are using Python 3 you will also need to modify the line, 
```python
PyTransSetup.compileName("QuartAxNC",True)
``` 
to, 
```python
PyTransSetup.compileName3("QuartAxNC",True)
```
in the setup file.

Back in shell, run the script,
```sh
python QuartAxNCsetup.py
```
which will install the `QuartAxNC` PyTranport library which can later be imported into a python script, e.g. `SimpleExample.py`.

# Contibutions, Pull requests and Issues
Any third party wishing to contribute to the PyTransport project is welcome to do so. We will attempt to implement submitted pull requests and fix issues submitted to this repository. If assistance or support is needed you can also contact us by email at j.ronayne@qmul.ac.uk  and d.mulryne@qmul.ac.uk.

# Licencing #
PyTransport is distributed under the GNU General Public License version 3, or (at your option) any later version. This license is bundled with the source code as LICENSE.txt.

Please visit the [PyTransport website](https://transportmethod.com) for further information and links to other repositories ultisiling the transport method.
