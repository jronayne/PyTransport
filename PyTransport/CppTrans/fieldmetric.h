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
nP=3;
	
   }
	
	
	//calculates fieldmetic()
	vector<double> fmetric(vector<double> f, vector<double> p)
	{
		vector<double> FM((2*nF)*(2*nF),0.0) ;
        
// metric
  double x0 = std::pow(p[2], 2.0);
  double x1 = 1.0/x0;
  double x2 = std::pow(std::sin(f[0]), 2.0);

 FM[0]=x1;

 FM[1]=0;

 FM[2]=1;

 FM[3]=0;

 FM[4]=0;

 FM[5]=x1/x2;

 FM[6]=0;

 FM[7]=1;

 FM[8]=1;

 FM[9]=0;

 FM[10]=x0;

 FM[11]=0;

 FM[12]=0;

 FM[13]=1;

 FM[14]=0;

 FM[15]=x0*x2;

         return FM;
	}
	
	
	
	//calculates ChristoffelSymbole()
	vector<double> Chroff(vector<double> f, vector<double> p)
	{
		vector<double> CS((2*nF)*(2*nF)*(2*nF),0.0);
	
// Christoffel
  double x0 = std::pow(std::sin(f[0]), 1.0);
  double x1 = 1.0*std::cos(f[0]);
  double x2 = x1/x0;

 CS[0]=0;

 CS[1]=0;

 CS[2]=0;

 CS[3]=0;

 CS[4]=0;

 CS[5]=0;

 CS[6]=0;

 CS[7]=0;

 CS[8]=0;

 CS[9]=0;

 CS[10]=0;

 CS[11]=0;

 CS[12]=0;

 CS[13]=0;

 CS[14]=0;

 CS[15]=-x0*x1;

 CS[16]=0;

 CS[17]=0;

 CS[18]=0;

 CS[19]=0;

 CS[20]=0;

 CS[21]=0;

 CS[22]=0;

 CS[23]=0;

 CS[24]=0;

 CS[25]=0;

 CS[26]=0;

 CS[27]=x2;

 CS[28]=0;

 CS[29]=0;

 CS[30]=x2;

 CS[31]=0;

 CS[32]=0;

 CS[33]=0;

 CS[34]=0;

 CS[35]=0;

 CS[36]=0;

 CS[37]=0;

 CS[38]=0;

 CS[39]=0;

 CS[40]=0;

 CS[41]=0;

 CS[42]=0;

 CS[43]=0;

 CS[44]=0;

 CS[45]=0;

 CS[46]=0;

 CS[47]=0;

 CS[48]=0;

 CS[49]=0;

 CS[50]=0;

 CS[51]=0;

 CS[52]=0;

 CS[53]=0;

 CS[54]=0;

 CS[55]=0;

 CS[56]=0;

 CS[57]=0;

 CS[58]=0;

 CS[59]=0;

 CS[60]=0;

 CS[61]=0;

 CS[62]=0;

 CS[63]=0;
        
		return CS;
	}
    

	
	// calculates RiemannTensor()
	vector<double> Riemn(vector<double> f, vector<double> p)
	{
		vector<double> RM((nF)*(nF)*(nF)*(nF),0.0);
		
// Riemann
  double x0 = 1.0*std::pow(p[2], 2.0)*std::pow(std::sin(f[0]), 2.0);
  double x1 = -x0;

 RM[0]=0;

 RM[1]=0;

 RM[2]=0;

 RM[3]=0;

 RM[4]=0;

 RM[5]=x0;

 RM[6]=x1;

 RM[7]=0;

 RM[8]=0;

 RM[9]=x1;

 RM[10]=x0;

 RM[11]=0;

 RM[12]=0;

 RM[13]=0;

 RM[14]=0;

 RM[15]=0;
     
        return RM;
	}

	// calculates RiemannTensor() covariant derivatives
	vector<double> Riemncd(vector<double> f, vector<double> p)
	{
		vector<double> RMcd((nF)*(nF)*(nF)*(nF)*(nF),0.0);
		
// Riemanncd

 RMcd[0]=0;

 RMcd[1]=0;

 RMcd[2]=0;

 RMcd[3]=0;

 RMcd[4]=0;

 RMcd[5]=0;

 RMcd[6]=0;

 RMcd[7]=0;

 RMcd[8]=0;

 RMcd[9]=0;

 RMcd[10]=0;

 RMcd[11]=0;

 RMcd[12]=0;

 RMcd[13]=0;

 RMcd[14]=0;

 RMcd[15]=0;

 RMcd[16]=0;

 RMcd[17]=0;

 RMcd[18]=0;

 RMcd[19]=0;

 RMcd[20]=0;

 RMcd[21]=0;

 RMcd[22]=0;

 RMcd[23]=0;

 RMcd[24]=0;

 RMcd[25]=0;

 RMcd[26]=0;

 RMcd[27]=0;

 RMcd[28]=0;

 RMcd[29]=0;

 RMcd[30]=0;

 RMcd[31]=0;
     
        return RMcd;
	}
    
    int getnF()
    {
        return nF;
    }
    


};
#endif

