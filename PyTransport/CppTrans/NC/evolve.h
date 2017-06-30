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

#include "fieldmetric.h"

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
vector<double> Dcovsig(vector<double> f,vector<double> p, double N)
	{
		model m;
		int nF = m.getnF();
		int nP = m.getnP();
		vector<double> covsigout(2*nF*2*nF);
		fieldmetric fmet;
		double Hi;
		double s=m.scale(f,p,N);
		double ds=m.dscale(f,p,N);
		Hi=m.H(f,p);

		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		vector<double> CHR;
		CHR = fmet.Chroff(f,p);
		double sum1=0.0;
		double sum2=0.0;
		double H=Hi;
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++)
        {	
			sum1=0.0;
			sum2=0.0;
            for(int m=0;m<nF;m++){
            for(int l=0;l<nF;l++){
            for(int xx=0;xx<nF;xx++){
				//sum2= sum2 - FMi[(2*nF)*(i)+(l)]*FMi[(2*nF)*(j+nF)+(xx+nF)]*CHR[(2*nF)*(2*nF)*(xx)+(2*nF)*(l+nF)+m+nF]*f[nF+m];
			}
			}
				sum1= sum1 + CHR[(2*nF)*(2*nF)*(i)+(2*nF)*(j+nF)+m+nF]*f[nF+m];
			}
			
			covsigout[i+ j*2*nF]=sum1/Hi;//+1./2.*sum2/Hi;
			covsigout[i+nF+(j)*2*nF]=0.0;
			covsigout[i+(j+nF)*2*nF]=0.0;
			covsigout[i+nF+(j+nF)*2*nF]=sum1/Hi;//+1./2.*sum2/Hi;
		}
		}
		return covsigout;
	}	

void evolveSig( double N, double yin[], double yp[], double paramsIn[])
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
	double Hi;
	Hi=m.H(fields,p);

	vector<double> dcov;
	dcov=Dcovsig(fields,p,N);
    
	
    
	for(int i=0;i<2*nF;i++){yp[i] = u1[i];}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
			double sum1=0.0;
            for(int m=0;m<2*nF;m++)
            {
				//if(i>k)
				sum1 = sum1  + 1./2.*dcov[i+ (m)*2*nF]*yin[2*nF+j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF+m+2*nF*i] + 1./2.*dcov[i+ (m)*2*nF]*yin[2*nF+m+2*nF*j] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF+i+2*nF*m];
                sum = sum + u2[i+m*2*nF]*yin[2*nF+m+2*nF*j] + u2[j+m*2*nF]*yin[2*nF+i+2*nF*m];
            }
            yp[2*nF+i+2*nF*j]=sum - sum1;

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
	vector<double> CHR;
	fieldmetric fmet;
	vector<double> FMi;
	FMi = fmet.fmetric(fields,p);
	CHR = fmet.Chroff(fields,p);
    vector<double> dcov;
    dcov=Dcovsig(fields,p,N);
	double Hi;
	Hi=m.H(fields,p);
    u1=m.u(fields,p);
    u2a=m.u(fields,p,k1,N);
    u2b=m.u(fields,p,k2,N);
    u2c=m.u(fields,p,k3,N);
    u3a=m.u(fields,p,k1, k2, k3, N);
    u3b=m.u(fields,p,k2, k1, k3, N);
    u3c=m.u(fields,p,k3, k1, k2, N);
	
	
    for(int i=0; i<2*nF; i++){yp[i] = u1[i];}
    //Calculates the evolution of the two point function for each k mode individually
    for(int i=0; i<2*nF; i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
            double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {

				sum3 = sum3 + 1./2.*dcov[i+ m*2*nF]*yin[2*nF+m+2*nF*j]+ 1./2.*dcov[i+ m*2*nF]*yin[2*nF+j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF+i+2*nF*m]+ 1./2.*dcov[j+(m)*2*nF]*yin[2*nF+m+2*nF*i];
                sum = sum + u2a[i+m*2*nF]*yin[2*nF+m+2*nF*j] + u2a[j+m*2*nF]*yin[2*nF+m+2*nF*i] ;
            }
            
            yp[2*nF+i+2*nF*j]=sum -sum3;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
			double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {

				sum3 = sum3 + 1./2.*dcov[i+ m*2*nF]*yin[2*nF + (2*nF*2*nF) + m+2*nF*j]+ 1./2.*dcov[i+ m*2*nF]*yin[2*nF + (2*nF*2*nF) + j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + (2*nF*2*nF) + m+2*nF*i]+ 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + (2*nF*2*nF) + i+2*nF*m];
                sum = sum + u2b[i+m*2*nF]*yin[2*nF + (2*nF*2*nF) + m+2*nF*j] + u2b[j+m*2*nF]*yin[2*nF + (2*nF*2*nF) + m+2*nF*i];
            }
            yp[2*nF + (2*nF*2*nF) + i+2*nF*j]=sum - sum3;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
			double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {

				sum3 = sum3 + 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 2*(2*nF*2*nF) +  m+2*nF*j]+ 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 2*(2*nF*2*nF) +  j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 2*(2*nF*2*nF) + m+2*nF*i]+ 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 2*(2*nF*2*nF) + i+2*nF*m];
                sum = sum + u2c[i+m*2*nF]*yin[2*nF + 2*(2*nF*2*nF) +  m+2*nF*j] + u2c[j+m*2*nF]*yin[2*nF + 2*(2*nF*2*nF) + m+2*nF*i];
            }
            yp[2*nF + 2*(2*nF*2*nF) +i+2*nF*j]=sum - sum3;
        }}
    
    
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
			double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {

				sum3 = sum3 + 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 3*(2*nF*2*nF) +  m+2*nF*j]+ 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 3*(2*nF*2*nF) +  j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 3*(2*nF*2*nF) + i+2*nF*m]+ 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 3*(2*nF*2*nF) + m+2*nF*i];
                sum = sum + u2a[i+m*2*nF]*yin[2*nF + 3*(2*nF*2*nF) +  m+2*nF*j] + u2a[j+m*2*nF]*yin[2*nF + 3*(2*nF*2*nF) + i+2*nF*m];
            }
            yp[2*nF + 3*(2*nF*2*nF) +i+2*nF*j]=sum -sum3;
        }}
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
			double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {

				sum3 = sum3 + 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 4*(2*nF*2*nF) +  m+2*nF*j]+ 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 4*(2*nF*2*nF) +  j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 4*(2*nF*2*nF) + i+2*nF*m]+ 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 4*(2*nF*2*nF) + m+2*nF*i];
                sum = sum + u2b[i+m*2*nF]*yin[2*nF + 4*(2*nF*2*nF) +  m+2*nF*j] + u2b[j+m*2*nF]*yin[2*nF + 4*(2*nF*2*nF) + i+2*nF*m];
            }
            yp[2*nF + 4*(2*nF*2*nF) +i+2*nF*j]=sum - sum3;
        }}
    
    
    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++)
        {
            double sum=0.0;
			double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {

                sum3 = sum3 + 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 5*(2*nF*2*nF) +  m+2*nF*j] + 1./2.*dcov[i+ m*2*nF]*yin[2*nF + 5*(2*nF*2*nF) +  j+2*nF*m] + 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 5*(2*nF*2*nF) + i+2*nF*m]+ 1./2.*dcov[j+(m)*2*nF]*yin[2*nF + 5*(2*nF*2*nF) + m+2*nF*i];
				sum = sum + u2c[i+m*2*nF]*yin[2*nF + 5*(2*nF*2*nF) +  m+2*nF*j] + u2c[j+m*2*nF]*yin[2*nF + 5*(2*nF*2*nF) + i+2*nF*m];
            }
            yp[2*nF + 5*(2*nF*2*nF) +i+2*nF*j]=sum - sum3;
        }}
    
	//Calculates the evolution of the three point function

    for(int i=0;i<2*nF;i++){for(int j=0;j<2*nF;j++){for(int k=0;k<2*nF;k++)
        {
            double sum=0.0;
            double sum2=0.0;
			double sum3=0.0;
            for(int m=0;m<2*nF;m++)
            {
				
				sum3 =  1./6.*dcov[j+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + k*2*nF + i*2*nF*2*nF]
							+ 1./6.*dcov[i+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + j*2*nF + k*2*nF*2*nF]
							+ 1./6.*dcov[i+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + k*2*nF + j*2*nF*2*nF]
							+ 1./6.*dcov[j+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + i*2*nF + k*2*nF*2*nF]
							+ 1./6.*dcov[k+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + i*2*nF + j*2*nF*2*nF]
							+ 1./6.*dcov[k+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + j*2*nF + i*2*nF*2*nF]
							+ 1./6.*dcov[j+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  i + m*2*nF + k*2*nF*2*nF]
							+ 1./6.*dcov[j+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  k + m*2*nF + i*2*nF*2*nF]
							+ 1./6.*dcov[k+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  j + m*2*nF + i*2*nF*2*nF]
							+ 1./6.*dcov[k+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  i + m*2*nF + j*2*nF*2*nF]
							+ 1./6.*dcov[i+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  j + m*2*nF + k*2*nF*2*nF]
							+ 1./6.*dcov[i+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  k + m*2*nF + j*2*nF*2*nF]
							+ 1./6.*dcov[k+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  i + j*2*nF + m*2*nF*2*nF]
							+ 1./6.*dcov[k+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  j + i*2*nF + m*2*nF*2*nF]							
							+ 1./6.*dcov[j+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  i + k*2*nF + m*2*nF*2*nF]
							+ 1./6.*dcov[j+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  k + i*2*nF + m*2*nF*2*nF]
							+ 1./6.*dcov[i+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  k + j*2*nF + m*2*nF*2*nF]
							+ 1./6.*dcov[i+ m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  j + k*2*nF + m*2*nF*2*nF];
							
                sum = sum - sum3 + u2a[i+m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)  +  m + j*2*nF + k*2*nF*2*nF] + u2b[j+m*2*nF]*yin[2*nF + 6*(2*nF*2*nF)
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


