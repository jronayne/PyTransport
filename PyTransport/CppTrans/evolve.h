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



//####################### File contains the evolution equations for the background and transport system in the correct form to be used by the integrateor ##########################

#include <iostream>
#include "moments.h"
#include "potential.h"
#include "model.h"
#include <math.h>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <time.h>

using namespace std;

//takes in current state of background system, y -- fields and field derivatives -- calculates dy/dN
void evolveB(double N, double yin[], double yp[], double paramsIn[])
{
    model m;
    int nF = m.getnF();
    int nP = m.getnP();
    vector<double> p(paramsIn,paramsIn+nP);
    vector<double> fields(yin, yin+2*nF);
    vector<double> ypin = m.u(fields,p);
    for (int i=0;i<2*nF;i++)
    {
       yp [i] = ypin[i];
    }
  //  delete [] in;
}

//takes in current state of background and 2pt transport system, y -- calculates dy/dN
void evolveSig( double N,  double yin[], double yp[], double paramsIn[])
{
    model m;
    double k;
    int nP = m.getnP();
    vector<double> p(paramsIn,paramsIn+nP);
    k=paramsIn[nP];
    int nF=m.getnF();
    vector<double> fields(yin, yin+2*nF);
    vector<double> u1=m.u(fields,p);
    vector<double> u2=m.u(fields,p,k,N);
    
    for(int i=0;i<2*nF;i++){yp[i] = u1[i];}
    
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2[i+m*2*nF]*yin[2*nF+m+2*nF*j] + u2[j+m*2*nF]*yin[2*nF+m+2*nF*i];
            }
            yp[2*nF+i+2*nF*j]=sum;
        }}
}


//takes in current state of background and 2pt Gamma transport system, y -- calculates dy/dN for
//models with field space metric

void evolveGam( double N, double yin[], double yp[], double paramsIn[])
{
    model m;
    double k;
	
    int nP = m.getnP();
	int nF = m.getnF();
    vector<double> p(paramsIn,paramsIn+nP);
    k=paramsIn[nP];
	int nT=1;
    vector<double> fields(yin, yin+2*nF);
    vector<double> u1=m.u(fields,p);
    vector<double> w2=m.w(fields,p,k,N);

    
	for(int i=0;i<2*nF;i++){yp[i] = u1[i];}
    
    for(int i=0;i<2*nT;i++){for(int j=0;j<2*nT;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nT;m++)
            {
				//if(i>k)
                sum = sum + w2[i+m*2*nT]*yin[2*nF+m+2*nT*j] + w2[j+m*2*nT]*yin[2*nF+i+2*nT*m];
            }
            yp[2*nF+i+2*nT*j]=sum ;

        }}

}

//takes in current state of background, the three 2pts functions needed to evolve the 3pt, and the 3pt, y, and calculates dy/dN
void evolveAlp(double N,  double yin[], double yp[], double paramsIn[])
{
    model m;
    double k1, k2, k3;
    int nF=m.getnF();
    int nP = m.getnP();
    vector<double> p(paramsIn,paramsIn+nP);
    k1=paramsIn[nP];
    k2=paramsIn[nP+1];
    k3=paramsIn[nP+2];
    vector<double> fields(yin, yin+2*nF);
    vector<double> u1, u2a, u2b, u2c, u3a, u3b, u3c;

    u1=m.u(fields,p);
    u2a=m.u(fields,p,k1,N);
    u2b=m.u(fields,p,k2,N);
    u2c=m.u(fields,p,k3,N);
    u3a=m.u(fields,p,k1, k2, k3, N);
    u3b=m.u(fields,p,k2, k1, k3, N);
    u3c=m.u(fields,p,k3, k1, k2, N);
    

    for(int i=0; i<2*nF; i++){yp[i] = u1[i];}
    
    for(int i=0; i<2*nF; i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2a[i+m*2*nF]*yin[2*nF+m+2*nF*j] + u2a[j+m*2*nF]*yin[2*nF+m+2*nF*i];
            }
            
            yp[2*nF+i+2*nF*j]=sum;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2b[i+m*2*nF]*yin[2*nF + (2*nF*2*nF) + m+2*nF*j] + u2b[j+m*2*nF]*yin[2*nF + (2*nF*2*nF) + m+2*nF*i];
            }
            yp[2*nF + (2*nF*2*nF) + i+2*nF*j]=sum;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2c[i+m*2*nF]*yin[2*nF + 2*(2*nF*2*nF) +  m+2*nF*j] + u2c[j+m*2*nF]*yin[2*nF + 2*(2*nF*2*nF) + m+2*nF*i];
            }
            yp[2*nF + 2*(2*nF*2*nF) +i+2*nF*j]=sum;
        }}
    
    
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2a[i+m*2*nF]*yin[2*nF + 3*(2*nF*2*nF) +  m+2*nF*j] + u2a[j+m*2*nF]*yin[2*nF + 3*(2*nF*2*nF) + i+2*nF*m];
            }
            yp[2*nF + 3*(2*nF*2*nF) +i+2*nF*j]=sum;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2b[i+m*2*nF]*yin[2*nF + 4*(2*nF*2*nF) +  m+2*nF*j] + u2b[j+m*2*nF]*yin[2*nF + 4*(2*nF*2*nF) + i+2*nF*m];
            }
            yp[2*nF + 4*(2*nF*2*nF) +i+2*nF*j]=sum;
        }}
    
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2c[i+m*2*nF]*yin[2*nF + 5*(2*nF*2*nF) +  m+2*nF*j] + u2c[j+m*2*nF]*yin[2*nF + 5*(2*nF*2*nF) + i+2*nF*m];
            }
            yp[2*nF + 5*(2*nF*2*nF) +i+2*nF*j]=sum;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++){for(int k=0;k<2*nF;k++)
        {
            double sum=0.0;
            double sum2=0.0;
            for(int m=0;m<2*nF;m++)
            {
                sum = sum + u2a[i+m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + j*2*nF + k*2*nF*2*nF] + u2b[j+m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)
                                                                                                                     +  i + m*2*nF + k*2*nF*2*nF] + u2c[k+m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  i + j*2*nF + m*2*nF*2*nF];
                for(int n=0;n<2*nF;n++){
                    sum2 = sum2 + u3a[i+ n*2*nF + m*2*nF*2*nF ]*yin[2*nF + 1*(2*nF*2*nF) + j + n*2*nF]*yin[2*nF +2* (2*nF*2*nF) + k + m*2*nF]
                    + u3b[j+ n*2*nF + m*2*nF*2*nF ]*yin[2*nF + 0*(2*nF*2*nF) + n + i*2*nF]*yin[2*nF +2* (2*nF*2*nF) + m + k*2*nF]
                    + u3c[k+ n*2*nF + m*2*nF*2*nF ]*yin[2*nF + 0*(2*nF*2*nF) + n+ i*2*nF]*yin[2*nF +1* (2*nF*2*nF) + j + m*2*nF]
                    - 1.*u3a[i+ n*2*nF + m*2*nF*2*nF ]*yin[2*nF + 4*(2*nF*2*nF) + n + j*2*nF]*yin[2*nF +5* (2*nF*2*nF) + m + k*2*nF]
                    - 1.*u3b[j+ n*2*nF + m*2*nF*2*nF ]*yin[2*nF + 3*(2*nF*2*nF) + i + n*2*nF]*yin[2*nF +5* (2*nF*2*nF) + m + k*2*nF]
                    - 1.*u3c[k+ n*2*nF + m*2*nF*2*nF ]*yin[2*nF + 3*(2*nF*2*nF) + i + n*2*nF]*yin[2*nF +4* (2*nF*2*nF) + j + m*2*nF];
                }}
            yp[2*nF + 6*(2*nF*2*nF)  +  i+2*nF*j+k*2*nF*2*nF]=sum+sum2;    

        }}}
    
    
}


