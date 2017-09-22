//#This file is part of PyTransport.

//#PyTransport is free software: you can redistribute it and/or modify
//#it under the terms of the GNU General Public License as published by
//#the Free Software Foundation, either version 3 of the License, or
//#(at your option) any later version.

//#PyTransport is distributed in the hope that it will be useful,
//#but WITHOUT ANY WARRANTY; without even the implied warranty of
//#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//#GNU General Public License for more details.

//#You should have received a copy of the GNU General Public License
//#along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.


// C++ file which defines the functions make available to Python through the MTeasy module.
#include <Python.h>
#include <iostream>
#include "numpy/arrayobject.h"

//don't adjust the labels at the end of the 4 lines below (they are used to fix directory structure)
#include"/home/jwr/Code/June/PyTransport2Dist/PyTransport/CppTrans/NC/evolve.h"//evolve
#include"/home/jwr/Code/June/PyTransport2Dist/PyTransport/CppTrans/NC/moments.h"//moments
#include"/home/jwr/Code/June/PyTransport2Dist/PyTransport/CppTrans/NC/model.h"//model
#include"/home/jwr/Code/June/PyTransport2Dist/PyTransport/CppTrans/stepper/rkf45.hpp"//stepper
//************************************************************************************************* 

#include <math.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <time.h>
# include <iomanip>
# include <cmath>

using namespace std;

// The line below is updated evey time the moduleSetup file is run.
// Package recompile attempted at: Tue Jun 20 15:41:14 2017


// Changes python array into C array (or rather points to pyarray data)
//    Assumes PyArray is contiguous in memory.             */
double *pyvector_to_Carray(PyArrayObject *arrayin)
{
    //  int i,n;
    //	n=arrayin->dimensions[0]; n is length of python array
    return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

int size_pyvector(PyArrayObject *arrayin)
{
      return arrayin->dimensions[0];  /* pointer to arrayin data as double */
}

// function to retun amplitude of potential
static PyObject* MT_V(PyObject* self,  PyObject *args)
{
    PyArrayObject *fieldsIn, *params;
    double *Cfields,*Cparams;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields = pyvector_to_Carray(fieldsIn);
    Cparams = pyvector_to_Carray(params);
    potential pp;
    int nF = pp.getnF(); if (nF!=size_pyvector(fieldsIn)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;}
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
        
    vector<double> vectIn;
    vectIn = vector<double>(Cfields, Cfields + nF);
    
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    return Py_BuildValue("d", pp.V(vectIn,Vparams));
}

// function to calculate derivatives of potential
static PyObject* MT_dV(PyObject* self,  PyObject *args)
{
    PyArrayObject *fieldsIn, *dVI, *params;
    double *Cfields, *dVC, *Cparams ;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields = pyvector_to_Carray(fieldsIn);
    Cparams = pyvector_to_Carray(params);
    
    potential pp;
    int nF = pp.getnF();if (nF!=size_pyvector(fieldsIn)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;}
    vector<double> vectIn;
    npy_intp dims[1];
    dims[0]=nF;
    
    dVI = (PyArrayObject*) PyArray_SimpleNew(1,dims,NPY_DOUBLE);
    dVC = (double*) dVI->data;
    
    vectIn = vector<double>(Cfields, Cfields +  nF);
    
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    vector<double> dVect = pp.dV(vectIn,Vparams);
    for(int i=0; i<nF;i++){dVC[i] = dVect[i];}
    
    return PyArray_Return(dVI);
}

// function to calculate derivatives of potential
static PyObject* MT_ddV(PyObject* self,  PyObject *args)
{
    PyArrayObject *fieldsIn, *ddVI, *params;
    double *Cfields, *ddVC, *Cparams ;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields = pyvector_to_Carray(fieldsIn);
    Cparams = pyvector_to_Carray(params);
    
    potential pp;
    int nF = pp.getnF();if (nF!=size_pyvector(fieldsIn)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;}
    
    vector<double> vectIn;
    npy_intp dims[2];
    dims[0]=nF; dims[1]=nF;
    
    ddVI = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    ddVC = (double*) ddVI->data;
    
    vectIn = vector<double>(Cfields, Cfields +  nF);
    
    int nP = pp.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    vector<double> ddVect = pp.dVV(vectIn,Vparams);
    for(int i=0; i<nF;i++){for(int j=0; j<nF;j++){ddVC[i+j*nF] = ddVect[i+j*nF];}}
   
    return PyArray_Return(ddVI);
}

// function to calculate Hubble rate
static PyObject* MT_H(PyObject* self,  PyObject *args)
{
    PyArrayObject *fields_dfieldsIn, *params;
    double *Cfields_dfields, *Cparams;
    if (!PyArg_ParseTuple(args, "O!O!",  &PyArray_Type, &fields_dfieldsIn,&PyArray_Type,&params)) {
        return NULL;}
    Cfields_dfields = pyvector_to_Carray(fields_dfieldsIn);
    model mm;
    int nF = mm.getnF(); if (2*nF!=size_pyvector(fields_dfieldsIn)){cout<< "\n \n \n field space array not of correct length\n \n \n ";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(Cfields_dfields, Cfields_dfields + 2*nF);
    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    return Py_BuildValue("d", mm.H(vectIn, Vparams));
}

// function to return number of fields (useful for cross checks)
static PyObject* MT_fieldNumber(PyObject* self,  PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;}
    model mm;
    return Py_BuildValue("i",mm.getnF());
}

// function to return number of parameters (useful for cross checks)
static PyObject* MT_paramNumber(PyObject* self,  PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;}
    model mm;
    return Py_BuildValue("i",mm.getnP());
}


// detect end of inflation within a suitable search window
static PyObject* MT_findEndOfInflation(PyObject* self, PyObject* args)
{
    PyArrayObject* initialCs;
    PyArrayObject* params;
    PyArrayObject* tols;
    double Ninit = 0.0;         // value of N corresponding to initial conditions
    double DeltaN = 10000.0;    // number of e-folds to search through

    // parse arguments; final argument is optional
    if(!PyArg_ParseTuple(args, "O!O!O!d|d", &PyArray_Type, &initialCs, &PyArray_Type, &params, &PyArray_Type, &tols,
                         &Ninit, &DeltaN))
      return NULL;

    // convert requested tolerances to a C array and extract absolute & relative error targets
    double* tolsC = pyvector_to_Carray(tols);

    double abserr = 1E-8;
    double relerr = 1E-8;
    if(size_pyvector(tols) != 2)
      {
        cout << "\n \n \n incorrect tolerances input, using defaults \n \n \n";
      }
    else
      {
        abserr = tolsC[0];
        relerr = tolsC[1];
      }

    // convert initial conditions to a C array
    double* CinitialCs = pyvector_to_Carray(initialCs);

    // check whether the expected number of initial conditions were supplied
    model mm;
    int nF = mm.getnF();
    if(size_pyvector(initialCs) != 2*nF)
      {
        cout << "\n \n \n field space array not of correct length \n \n \n";
        Py_RETURN_NONE;
      }

    // convert parameter list to a C array
    double* Cparams = pyvector_to_Carray(params);
    int nP = mm.getnP();
    if(size_pyvector(params) != nP)
      {
        cout << "\n \n \n parameters array not of correct length \n \n \n";
        Py_RETURN_NONE;
      }

    // allocate working space for the stepper
    // TODO: consider absorbing these allocations within a janitor object
    double* y = new double[2*nF];       // current values
    double* dy = new double[2*nF];      // derivatives

    // initialize y using supplied initial conditions
    for(int i = 0; i < 2*nF; ++i)
      {
        y[i] = CinitialCs[i];
      }

    // populate dy for initial step
    double N = Ninit;
    const double Nstop = Ninit + DeltaN;
    evolveB(N, y, dy, Cparams);

    // integrated background until we encounter the end-of-inflation, or the end of the search window
    int flag = -1;      // '-1' puts the integrator into 'single-step' mode, so it returns after taking one stride
    while(N < Nstop)
      {
        flag = r8_rkf45(evolveB, 2*nF, y, dy, &N, Nstop, &relerr, abserr, flag, Cparams);

        // detect some error conditions
        if(flag == 50)
          {
            cout << "\n \n \n Integrator failed at time N = " <<N <<" \n \n \n";
            return Py_BuildValue("d", N);
          }

        // compute value of epsilon
        vector<double> vecy(y, y+2*nF);
        vector<double> vecParams(Cparams, Cparams+nP);

        double eps = mm.Ep(vecy, vecParams);

        // break out of search if we have an unacceptable result, or if we are now past the end of inflation
        if(eps < 0 || eps > 1)
          {
            // TODO: consider absorbing in a janitor
            delete[] y;
            delete[] dy;
            return Py_BuildValue("d", N);
          }

        flag = -2;      // '-2' means continue as normal in 'single-step' mode
      }

    // deallocate workspace
    // TODO: consider absorbing in a janitor
    delete[] y;
    delete[] dy;

    Py_RETURN_NONE;
}


// function to calculate background evolution
static PyObject* MT_backEvolve(PyObject* self,  PyObject *args)
{
    PyArrayObject *initialCs, *t, *backOut,*backOutT, *params, *tols;
    double *CinitialCs, *tc, *Cparams, *tolsC ;
    bool exit;
    
    if (!PyArg_ParseTuple(args, "O!O!O!O!b",&PyArray_Type, &t, &PyArray_Type, &initialCs, &PyArray_Type, &params, &PyArray_Type, &tols, &exit)) {
        return NULL;}
    
    tolsC = pyvector_to_Carray(tols);
    double abserr, relerr;
    if (2!=size_pyvector(tols)){cout<< "\n \n \n incorrect tolorances input, using defaults  \n \n \n";
        abserr = pow(10,-8.); relerr = pow(10,-8.);}
    else {
        abserr =tolsC[0];relerr = tolsC[1];}
    
    
    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    model mm;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "\n \n \n field space array not of correct length \n \n \n";    Py_RETURN_NONE;
    }


    
    double N=tc[0];
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs + 2*nF);
    back b(nF, vectIn);
    
       int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    
    
    int flag=-1;
    double *y; y = new double[2*nF];
    double *yp; yp= new double[2*nF];
    
    for (int i=0;i<2*nF;i++){y[i] = CinitialCs[i];}
    
    
    if (exit == false){
    int nt = t->dimensions[0];
    
    npy_intp dims[2];
    dims[1]=1+2*nF; dims[0]=nt;
    double * backOutC;
    backOut = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    backOutC = (double *) PyArray_DATA(backOut);
    evolveB(N, y, yp, Cparams);
    // run background *********************
    for(int ii=0; ii<nt; ii++ ){
        while (N<tc[ii]){
            flag = r8_rkf45(evolveB , 2*nF, y, yp, &N, tc[ii], &relerr, abserr, flag, Cparams );
            if (flag== 50){cout<< "\n \n \n Integrator failed at time N = " <<N <<" \n \n \n";  return Py_BuildValue("d", N);}
            flag=-2;
        }
        backOutC[ii*(2*nF+1)]=N;
        for(int i=0;i< 2*nF;i++){
            backOutC[ii*(2*nF+1)+i+1]=y[i];} // output array
    }
    }
    
    if (exit == true){
        int nt = t->dimensions[0];
        vector<double> vecy;
        vector<double> Vparams;
        npy_intp dims[2];
        dims[1]=1+2*nF; dims[0]=nt;
        double * backOutCT;
        backOutT = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
        backOutCT = (double *) PyArray_DATA(backOutT);
        evolveB(N, y, yp, Cparams);
        // run background *********************
        {int ii =0;double eps=0.0;
        while (eps<1 && ii<nt){
            while (N<tc[ii]){
                flag = r8_rkf45(evolveB , 2*nF, y, yp, &N, tc[ii], &relerr, abserr, flag, Cparams );
                if (flag== 50){cout<< "\n \n \n Integrator failed at time N = " <<N <<" \n \n \n"; return Py_BuildValue("d", N);}
                flag=-2;
            }
            backOutCT[ii*(2*nF+1)]=N;
            for(int i=0;i< 2*nF;i++){
                backOutCT[ii*(2*nF+1)+i+1]=y[i];} // outputs to file at each step
            
        vecy = vector<double>(y, y + 2*nF);
        Vparams = vector<double>(Cparams, Cparams +  nP);
            eps = mm.Ep(vecy,Vparams);
        ii = ii+1;
        }
            cout << ii <<endl;
        npy_intp dims2[2];
        dims2[1]=1+2*nF; dims2[0]=ii;
        double * backOutC;
        backOut = (PyArrayObject*) PyArray_SimpleNew(2,dims2,NPY_DOUBLE);
        backOutC = (double *) PyArray_DATA(backOut);
     
        for(int jj = 0; jj<ii; jj++){backOutC[jj*(2*nF+1)]=tc[jj];
        for(int i=0;i< 2*nF;i++){
            backOutC[jj*(2*nF+1)+i+1]=backOutCT[jj*(2*nF+1)+i+1] ;}}
        }
        }
    
    
    
    delete[] y; delete[] yp;
    return PyArray_Return(backOut);
}

// function to calculate 2pt evolution
static PyObject* MT_sigEvolve(PyObject* self,  PyObject *args)
{
    PyArrayObject *initialCs, *t, *sigOut, *params, *tols;
    double *CinitialCs, *tc, k, *Cparams, *tolsC;
    bool full;
    if (!PyArg_ParseTuple(args, "O!dO!O!O!b", &PyArray_Type, &t, &k, &PyArray_Type, &initialCs,&PyArray_Type, &params, &PyArray_Type, &tols,&full)) {
        return NULL;}
    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    
//    if (full != 0 && full !=1 ){ full=1; cout << "\n \n \n Number out of range, defaulted to full outout mode \n \n \n";}
    
    tolsC = pyvector_to_Carray(tols);
    double rtol, atol;
    if (2!=size_pyvector(tols)){cout<< "\n \n \n incorrect tolorances input, using defaults  \n \n \n";
        atol = pow(10,-8.); rtol = pow(10,-8.);}
    else {
        atol =tolsC[0];rtol = tolsC[1];}


    model mm;
    potential pott;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "\n \n \n field space array not of correct length, not proceeding further \n \n \n";    Py_RETURN_NONE;}
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs + 2*nF);
    
    int nP = mm.getnP(); if (nP!=size_pyvector(params)){cout<< "\n \n \n parameters array not of correct length, not proceeding further \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
    
    // we use a scaling below that we rescale back at the end (so the final answer is as if the scaling was never there -- this helps standarise the rtol and atol needed for the same model run with differnet initial conditions
    double kn = 1.0; 
    double kscale = k;    
    double Nstart=tc[0] - log(kscale);
    
    
    sigma sig(nF, kn, Nstart, vectIn, Vparams) ; // instance of sigma object which fixs ics
    
    double* y; // set up array for ics
    y = new double[2*nF + 2*nF*2*nF];
    
    for(int i=0; i<2*nF;i++){y[i] = CinitialCs[i];} // fix values of input array
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF+ i+2*nF*j] = sig.getS(i,j);}}
    
    double* paramsIn; // array of parameters to pass to LHS of ODE routine
    
    paramsIn = new double[1+nP];
    for(int i=0; i<nP;i++) paramsIn[i]=Vparams[i];
    paramsIn[nP]=kn;
    
    // evolve a 2pt run **************************
    
    double N=Nstart;
    double* yp ; yp = new double [2*nF + 2*nF*2*nF];
    vector<double> Ni;
    double zz=0;
    
    int flag=-1;
    evolveSig(N, y, yp, paramsIn);
    vector<double> fieldIn(2*nF);
    fieldIn = vector<double>(y,y+2*nF);
    Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
    zz=0;
    for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
        zz=zz+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];}
    }
    
    int nt = t->dimensions[0];
    
    int size;
    if (full ==true){size = 1+2*nF + 1+ 2*nF*2*nF;}
    if (full ==false){size = 1 + 1;}
    
    npy_intp dims[2];
    dims[1]=size; dims[0]=nt;
    double * sigOutC;
    sigOut = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    sigOutC = (double *) PyArray_DATA(sigOut);
    
    
    for(int ii=0; ii<nt; ii++ ){
        while (N<tc[ii]-log(kscale)){
            flag = r8_rkf45(evolveSig , 2*nF+2*nF*2*nF, y, yp, &N, tc[ii]-log(kscale), &rtol, atol, flag, paramsIn );
            if (flag== 50){cout<< "\n \n \n Integrator failed at time N = " <<N <<" \n \n \n"; return Py_BuildValue("d", N);}
            flag = -2;
        }
        fieldIn = vector<double>(y,y+2*nF);
        
        sigOutC[ii*size] = N+log(kscale);
        
        Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
        zz=0;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
            zz=zz+Ni[i]*Ni[j]*y[2*nF + i + j*2*nF];}}
        
        sigOutC[ii*size+1] = zz/kscale/kscale/kscale;
        
        
        if(full==true){
            for(int i=0;i<2*nF;i++)
            {
                sigOutC[ii*(size)+i+2]=y[i];
            }
            for(int i=2*nF;i<2*nF+ 2*nF*2*nF;i++)
            {
                sigOutC[ii*(size)+i+2]=y[i]/kscale/kscale/kscale;
            }
        }
        
    }
    
    delete [] y; delete [] yp;
    delete [] paramsIn;
    return PyArray_Return(sigOut);
}

// function to calculate 3pt evolution
static PyObject* MT_alphaEvolve(PyObject* self,  PyObject *args)
{
    PyArrayObject *initialCs, *t, *alpOut,*params, *tols;
    double k1, k2, k3, Nstart, *CinitialCs, *tc,*Cparams, *tolsC;
    bool full;
    if (!PyArg_ParseTuple(args, "O!dddO!O!O!b", &PyArray_Type, &t, &k1,&k2,&k3, &PyArray_Type, &initialCs,&PyArray_Type,&params,&PyArray_Type,&tols, &full)) {
        return NULL; }
    CinitialCs = pyvector_to_Carray(initialCs);
    tc = pyvector_to_Carray(t);
    int nt = t->dimensions[0];
    
    tolsC = pyvector_to_Carray(tols);
    double rtol, atol;
    if (2!=size_pyvector(tols)){cout<< "\n \n \n incorrect tolorances input, using defaults  \n \n \n";
        atol = pow(10,-8.); rtol = pow(10,-8.);}
    else {
        atol =tolsC[0];rtol = tolsC[1];}
 
 //   if (full != 0 && full !=1 ){ full=1;cout << "\n \n \n Number out of range, defaulted to full outout mode \n \n \n";}
    model mm;
    int nF=mm.getnF(); if (2*nF!=size_pyvector(initialCs)){cout<< "\n \n \n field space array not of correct length, not proceeding further \n \n \n";    Py_RETURN_NONE;}
    
    int nP = mm.getnP();if (nP!=size_pyvector(params)){cout<< "\n \n \n  parameters array not of correct length, not proceeding further \n \n \n";  Py_RETURN_NONE;}
    Cparams = pyvector_to_Carray(params);
    vector<double> Vparams; Vparams = vector<double>(Cparams, Cparams +  nP);
 
// we use a scaling below that we rescale back at the end (so the final answer is as if the scaling was never there -- this helps standarises the rtol and atol needed for the same model run with differnet initial conditions
    
    double kscale = (k1+k2+k3)/3.;
    double k1n = k1/kscale; double k2n = k2/kscale; double k3n = k3/kscale;
    Nstart=tc[0] -log(kscale);
    double N=Nstart; // reset N
    
    // do not alter the comment at the end of the next line -- used by preprocessor
    // ****************************************************************************
        
    vector<double> vectIn;
    vectIn = vector<double>(CinitialCs, CinitialCs+2*nF);
    
    sigma sig1(nF, k1n, Nstart, vectIn,Vparams)  ; // 3 instances of sigmmas
    sigma sig2(nF, k2n, Nstart, vectIn,Vparams)  ;
    sigma sig3(nF, k3n, Nstart, vectIn,Vparams)  ;
    sigmaI sig1I(nF, k1n, Nstart, vectIn,Vparams)  ; // 3 instances of sigma imaginary
    sigmaI sig2I(nF, k2n, Nstart, vectIn,Vparams)  ;
    sigmaI sig3I(nF, k3n, Nstart, vectIn,Vparams)  ;
    alpha alp(nF, k1n, k2n, k3n, Nstart, vectIn, Vparams); // instance of alpha
    
    double* y; // array for initial conditions
    y = new double[2*nF + 6*2*nF*2*nF + 2*nF*2*nF*2*nF];
    
    for(int i=0; i<2*nF;i++){y[i] = CinitialCs[i];}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF+ i+2*nF*j] = sig1.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 1*(2*nF*2*nF)+ i+2*nF*j] = sig2.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 2*(2*nF*2*nF)+ i+2*nF*j] = sig3.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 3*(2*nF*2*nF)+ i+2*nF*j] = sig1I.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 4*(2*nF*2*nF)+ i+2*nF*j] = sig2I.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){y[2*nF + 5*(2*nF*2*nF)+ i+2*nF*j] = sig3I.getS(i,j);}}
    for(int i=0; i< 2*nF;i++){for(int j=0;j<2*nF;j++){for(int k=0; k<2*nF;k++){y[2*nF + 6*(2*nF*2*nF)+ i+2*nF*j + 2*nF*2*nF*k] = alp.getA(i,j,k);}}}
    
    
    double* paramsIn2; // array for parameters of RHS of ODE routine
    paramsIn2 = new double[3+nP];
    for(int i=0; i<nP;i++) paramsIn2[i]=Vparams[i];
    paramsIn2[nP]=k1n;
    paramsIn2[nP+1]=k2n;
    paramsIn2[nP+2]=k3n;
    
    
    double ZZZ=0., ZZ1=0., ZZ2=0., ZZ3=0.; //  for zeta zeta calcs
    vector<double> Ni, Nii1, Nii2, Nii3 ; // for N transforms to get to zeta
    double *yp; yp=new double[2*nF +6*2*nF*2*nF+  2*nF*2*nF*2*nF];
    
    npy_intp dims[2];
    int size;
    if (full==false){size =   5;}
    if (full==true){size =  5+  2*nF + 6*2*nF*2*nF+2*nF*2*nF*2*nF;}
    dims[1]=size; dims[0]=nt;
    double * alpOutC;
    alpOut = (PyArrayObject*) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    alpOutC = (double *) PyArray_DATA(alpOut);
    
    
    evolveAlp(N, y, yp, paramsIn2);
    int flag=-1;
    
    // run alpha *******************************************
    vector<double> fieldIn(2*nF);
    
    for(int ii=0; ii<nt; ii++ ){
        while (N<(tc[ii]-log(kscale))){
            flag = r8_rkf45(evolveAlp, 2*nF + 6*(2*nF*2*nF) + 2*nF*2*nF*2*nF, y, yp, &N, tc[ii]-log(kscale), &rtol, atol, flag, paramsIn2);
            if (flag== 50){cout<< "\n \n \n Integrator failed at time N = " <<N <<" \n \n \n"; return Py_BuildValue("d", N);}
            flag=-2;
        }
        fieldIn = vector<double>(y,y+2*nF);
        Ni=mm.N1(fieldIn,Vparams,N); // calculate N,i array
        Nii1=mm.N2(fieldIn,Vparams,k1n,k2n,k3n,N); // claculate N,ij array for first arrangement of ks
        Nii2=mm.N2(fieldIn,Vparams,k2n,k1n,k3n,N); // for second
        Nii3=mm.N2(fieldIn,Vparams,k3n,k1n,k2n,N); // etc
        
        ZZ1=0.;
        ZZ2=0.;
        ZZ3=0.;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF; j++){
            ZZ1=ZZ1+Ni[i]*Ni[j]*(y[2*nF + i + j*2*nF] );
            ZZ2=ZZ2+Ni[i]*Ni[j]*y[2*nF + (2*nF*2*nF) + i + j*2*nF];
            ZZ3=ZZ3+Ni[i]*Ni[j]*y[2*nF + 2*(2*nF*2*nF) + i + j*2*nF];
        }}   


        ZZZ=0.;
        for(int i=0; i<2*nF;i++){for(int j=0; j<2*nF;j++){for(int k=0; k<2*nF;k++){
            ZZZ=ZZZ + Ni[i]*Ni[j]*Ni[k]*y[2*nF + 6*(2*nF*2*nF) + i + j*2*nF+ k*2*nF*2*nF];
            for(int l=0; l<2*nF;l++){ZZZ=ZZZ+(Nii1[i+j*2*nF]*Ni[k]*Ni[l]*y[2*nF + 1*(2*nF*2*nF) + i+k*2*nF]*y[2*nF+2*(2*nF*2*nF)+j+l*2*nF]
                                              +Nii2[i+j*2*nF]*Ni[k]*Ni[l]*y[2*nF + 0*(2*nF*2*nF) + i+k*2*nF]*y[2*nF + 2*(2*nF*2*nF) + j+l*2*nF]
                                              +Nii3[i+j*2*nF]*Ni[k]*Ni[l]*y[2*nF + 0*(2*nF*2*nF) + i+k*2*nF]*y[2*nF + 1*(2*nF*2*nF) + j+l*2*nF]);
            }}}}
        
        alpOutC[ii*size] =  N+log(kscale);
        //cout << N+log(kscale) << endl;
        alpOutC[ii*size+1] = ZZ1/kscale/kscale/kscale;
        alpOutC[ii*size+2] = ZZ2/kscale/kscale/kscale;
        alpOutC[ii*size+3] = ZZ3/kscale/kscale/kscale;
        alpOutC[ii*size+4] = ZZZ/kscale/kscale/kscale/kscale/kscale/kscale;
        
        if(full==true){
            for(int i=0;i<2*nF ;i++){
                alpOutC[ii*size+5+i] =  y[i] ;   }
            
            for(int i=2*nF;i<2*nF + 6*(2*nF*2*nF);i++){
                alpOutC[ii*size+5+i] =  y[i]/kscale/kscale/kscale ;   }

            for(int i=2*nF + 6*(2*nF*2*nF);i<2*nF + 6*(2*nF*2*nF)+ 2*nF*2*nF*2*nF;i++){
                alpOutC[ii*size+5+i] =  y[i]/kscale/kscale/kscale/kscale/kscale/kscale ;   }
        }
        
    }
    
    delete [] y;  delete [] paramsIn2; delete [] yp;
    
    return PyArray_Return(alpOut);
}


static char PyTrans_docs[] =
"This is PyTrans, a package for solving the moment transport equations of inflationary cosmology\n";

// **************************************************************************************
static PyMethodDef PyTransDQuadNC_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs}, {"findEndOfInflation", (PyCFunction)MT_findEndOfInflation,        METH_VARARGS, PyTrans_docs}, {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},    {NULL}};//FuncDef
// do not alter the comment at the end of preceeding line -- it is used by preprocessor

#ifdef __cplusplus
extern "C" {
#endif

// **************************************************************************************    
     //modDef
// do not alter the comment at the end of preceeding line -- it is used by preprocessor
    
// **************************************************************************************
void initPyTransDQuadNC(void)    {        Py_InitModule3("PyTransDQuadNC", PyTransDQuadNC_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc
// do not alter the comment at the end of preceeding line -- it is used by preprocessor

#ifdef __cplusplus
}
#endif
