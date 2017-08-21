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

// This file contains a prototype of the potential.h file of PyTransport -- it is edited by the PyTransScripts module

#ifndef FIELDMETRIC_H  // Prevents the class being re-defined
#define FIELDMETRIC_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;


class fieldmetric
{
private:
	int nF; // field number
    int nP; // params number which definFs potential

public:
    fieldmetric()
   {
// #FP
nF=2;
nP=4;
	
   }
	
	
	//calculates fieldmetic()
	vector<double> fmetric(vector<double> f, vector<double> p)
	{
		vector<double> FM((2*nF)*(2*nF),0.0) ;
        
// metric

 FM[0]=1.0*pow(p[3], -2.0);

 FM[2]=1.00000000000000;

 FM[5]=1.0*pow(p[3], -2.0)*pow(sin(f[0]), -2.0);

 FM[7]=1.00000000000000;

 FM[8]=1.00000000000000;

 FM[10]=1.0*pow(p[3], 2.0);

 FM[13]=1.00000000000000;

 FM[15]=1.0*pow(p[3], 2.0)*pow(sin(f[0]), 2.0);

         return FM;
	}
	
	
	
	//calculates ChristoffelSymbole()
	vector<double> Chroff(vector<double> f, vector<double> p)
	{
		vector<double> CS((2*nF)*(2*nF)*(2*nF),0.0);
	
// Christoffel

 CS[15]=-1.0*pow(sin(f[0]), 1.0)*cos(f[0]);

 CS[27]=1.0*1.0/sin(f[0])*cos(f[0]);

 CS[30]=1.0*1.0/sin(f[0])*cos(f[0]);
        
		return CS;
	}
    

	
	// calculates RiemannTensor()
	vector<double> Riemn(vector<double> f, vector<double> p)
	{
		vector<double> RM((nF)*(nF)*(nF)*(nF),0.0);
		
// Riemann

 RM[5]=1.0*pow(p[3], 2.0)*pow(sin(f[0]), 2.0);

 RM[6]=-1.0*pow(p[3], 2.0)*pow(sin(f[0]), 2.0);

 RM[9]=-1.0*pow(p[3], 2.0)*pow(sin(f[0]), 2.0);

 RM[10]=1.0*pow(p[3], 2.0)*pow(sin(f[0]), 2.0);
     
        return RM;
	}

	// calculates RiemannTensor() covariant derivatives
	vector<double> Riemncd(vector<double> f, vector<double> p)
	{
		vector<double> RMcd((nF)*(nF)*(nF)*(nF)*(nF),0.0);
		
// Riemanncd
     
        return RMcd;
	}
    
    int getnF()
    {
        return nF;
    }
    


};
#endif

