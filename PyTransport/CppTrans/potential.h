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

#ifndef POTENTIAL_H  // Prevents the class being re-defined
#define POTENTIAL_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;

// #Rewrite
// Potential file rewriten at Wed Apr 24 14:30:20 2019

class potential
{
private:
	int nF; // field number
	int nP; // params number which definFs potential
    
    
public:
	// flow constructor
	potential()
	{
// #FP
nF=2;
nP=3;

//        p.resize(nP);
        
// pdef

    }
	
    //void setP(vector<double> pin){
    //    p=pin;
    //}
	//calculates V()
	double V(vector<double> f, vector<double> p)
	{
		double sum ;
        
// Pot
  sum=0.5*std::pow(f[0], 2.0)*std::pow(p[0], 2.0) + 0.5*std::pow(f[1], 2.0)*std::pow(p[1], 2.0);
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);
	
// dPot

 sum[0]=1.0*std::pow(f[0], 1.0)*std::pow(p[0], 2.0);

 sum[1]=1.0*std::pow(f[1], 1.0)*std::pow(p[1], 2.0);
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot
  double x0 = 1.0*std::pow(p[0], 2.0);
  double x1 = std::pow(std::sin(f[0]), 1.0);
  double x2 = std::cos(f[0]);
  double x3 = 1.0*std::pow(p[1], 2.0);
  double x4 = -std::pow(f[1], 1.0)*x2*x3/x1;

 sum[0]=x0;

 sum[2]=x4;

 sum[1]=x4;

 sum[3]=std::pow(f[0], 1.0)*x0*x1*x2 + x3;
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot
  double x0 = std::pow(f[1], 1.0);
  double x1 = std::sin(f[0]);
  double x2 = std::pow(x1, 2.0);
  double x3 = 1.0/x2;
  double x4 = std::cos(f[0]);
  double x5 = std::pow(x4, 2);
  double x6 = 1.0*x5;
  double x7 = x3*x6;
  double x8 = std::pow(p[1], 2.0);
  double x9 = 1.0*x8;
  double x10 = x0*x8;
  double x11 = -x0*x9*(-x7 - 1.0) + x10*x7;
  double x12 = std::pow(x1, 1.0);
  double x13 = 1.0*std::pow(p[0], 2.0);
  double x14 = x12*x13*x4;
  double x15 = std::pow(f[0], 1.0);
  double x16 = x4/x12;
  double x17 = x16*(x14*x15 + x9);
  double x18 = 2.0*x10*x5;
  double x19 = x14 - x16*x9 - 1.0*x17;

 sum[0]=0;

 sum[4]=x18*x3;

 sum[2]=x11;

 sum[6]=x19;

 sum[1]=x11;

 sum[5]=x19;

 sum[3]=-x13*x15*(1.0*x2 - x6) + x14 - 2.0*x17;

 sum[7]=-x18;
       
        return sum;
	}
    
    int getnF()
    {
        return nF;
    }
    
    int getnP()
    {
        return nP;
    }

};
#endif