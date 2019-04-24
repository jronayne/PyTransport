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

//model class file contains the defining features of the model -- the u1, u2 flow tesors, and the A B C tesors and u3 flow tensor as well as the guage transform N tensors  
#ifndef MODEL_H  // Prevents the class being re-defined
#define MODEL_H 

#include "fieldmetric.h"
#include "potential.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;

class model
{
private:
	int nF;  // field number
	int nP; // params number which definFs potential
	potential pot; // potential which defines model
    fieldmetric fmet;

public:
	// constructor
	model()
	{
	   potential pot;
	   fieldmetric fmet;
        nP=pot.getnP();
        nF=pot.getnF();

    }
    
    // function returns Hubble rate
	double H(vector<double> f, vector<double> p   )
	{
		double Hi2;
		double Vi;
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		Vi=pot.V(f,p);
		Hi2=0.;

        
		for(int i=0; i<nF; i++)
		{	for(int j=0; j<nF; j++)
			{
				Hi2=Hi2+1./3.*(FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF+j]*f[nF+i]/2.);
			}
	
		}

		Hi2=Hi2 + 1./3.*Vi;
				
        return sqrt(Hi2);
	}

    // function returns H dot
    double Hdot(vector<double> f, vector<double> p)
	{
        double sum=0.;

		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		for(int i=0; i<nF; i++)
		{for(int j=0; j<nF; j++)
			{
				sum= sum - 1./2.*(FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF+i]*f[nF+j]);
				}
		
		}

		return sum;
	}
    // function returns a double dot
    double addot(vector<double> f, vector<double> p)
    {
        double sum=0.;
        double addot;
        double Vi;
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
        Vi=pot.V(f,p);
        for(int i=0; i<nF; i++)
        {for(int j=0; j<nF; j++)
			{
				sum= sum -1./2.*(FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF+i]*f[nF+j]);
			}
        }
        addot=-1./3.*(sum+Vi);
		
		
        return addot;
    }
    // function returns epsilon
	double Ep(vector<double> f,vector<double> p)
	{
		double Hi = H(f,p);
		double mdotH=0.;
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);

		for(int i=0; i<nF; i++)
		{for(int j=0; j<nF; j++)
			{
			mdotH= mdotH + 1./2.*(FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF+i]*f[nF+j]);// should this be a minus?
			}

		}
		
  		return mdotH/(Hi*Hi);
	}
    
    
    // function returns number of fields
    int getnF()
    {
        return nF;
    }
    
    // function returns number of fields
    int getnP()
    {
        return nP;
    }
	
	// a rescale function for delta dot phi to improve performance
    double scale(vector<double> f, vector<double> p,  double N)
    {
        double k=1.0;
        double a =exp(N);
        double Hi=H(f,p);
    //    return a;
    //    return a + (1.0-a)/pow(1. + k/a/Hi,2.0);
        return  a/(1.+a*Hi/k)/Hi; //+ 1./pow(1. + k/a/Hi,10.0);
    }
    
    // derivative of rescaling function divided by rescaling function
    double dscale(vector<double> f, vector<double> p, double N)
    {
        double k = 1.0;
        double a = exp(N);
        double Hi = H(f,p);
        double Hdi = Hdot(f,p);
        //return a*Hi;
       // return  a*Hi - a *Hi/pow(1. + k/a/Hi,2.0) + (1.0-a)*2.0*k/a/a/Hi/Hi*(a*Hi*Hi+a*Hdi)/pow(1. + k/a/Hi,3.0);
        return  -Hdi/Hi/Hi*a/(1.+a*Hi/k) + a/(1.+a*Hi/k) -a*(a*Hi*Hi/k + a*Hdi/k)/(1.+a*Hi/k)/(1.+a*Hi/k)/Hi; //+ 10.0*k/a/a/Hi/Hi*(a*Hi*Hi+a*Hdi)/pow(1. + k/a/Hi,11.0);
    }
    
    
    
    
    
		//ooooo  oooo  oo                      ooooooooooo                                                           
		// 888    88 o888                      88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		// 888    88  888       ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		// 888    88  888                          888    888          888   888          888 888     888 888        
		//  888oo88  o888o                        o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       

	
    
    // calculates u1 for models with field space metric (Christoffel terms are included here)
    
	vector<double> u(vector<double> f,vector<double> p)
	{
		vector<double> u1out(2*nF);
		vector<double> dVi;
		double Hi;
		Hi=H(f,p);
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);

		vector<double> CHR;
		CHR = fmet.Chroff(f,p);

		
		for(int i=0; i<nF; i++)
		{	
	
			u1out[i]  = f[nF+i]/Hi;
		}

		dVi=pot.dV(f,p);

		for(int i=0; i<nF; i++)			
		{	
			double sum =0.0;
			double sum1 =0.0;
			for(int j=0; j<nF; j++)
			{
				for(int k=0; k<nF; k++){
					sum1 =sum1 + CHR[(2*nF)*(2*nF)*i + (2*nF)*(j+nF)+k+nF]*f[nF+j]*f[nF+k];
				}
				sum = sum  +FMi[(2*nF)*(i)+(j)]*(-dVi[j])/Hi;
			}
			u1out[nF+i]  = sum -3.*Hi*f[nF+i]/Hi - sum1/Hi;

			}
		
		return u1out;
	}

    
	
		//ooooo  oooo  ooooooo                       ooooooooooo                                                           
		// 888    88 o88     888                     88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		// 888    88       o888       ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		// 888    88    o888   o                         888    888          888   888          888 888     888 888        
		//  888oo88  o8888oooo88                        o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       
	// calculates u2 tensor for models will non-trivial field space metric -- note that
    // Christoffel terms in the equations of motion are added in the evolve.h field

	vector<double> u(vector<double> f,vector<double> p, double k1, double N)
	{
		double a = exp(N);
		double ep = Ep(f,p);
		vector<double> u2out(2*nF*2*nF);
		double Hi=H(f,p);
        double s=scale(f,p,N);
        double ds=dscale(f,p,N);
		vector<double> dVVi;
		dVVi = pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		vector<double> RMi;
		RMi = fmet.Riemn(f,p);
		vector<double> u1 =u(f,p);
		
		
		double sum1 = 0.0;
		double sum2 = 0.0;
		double sum3 = 0.0;
		double sum4 = 0.0;

		
		for(int i = 0; i<nF; i++){for(int j = 0; j<nF; j++){
			for(int l=0; l<nF; l++)
			{for(int m=0; m<nF; m++)
			{for(int xx=0; xx<nF; xx++){
			sum1 = sum1 + 1./2.*FMi[(2*nF)*xx +i]*(RMi[(nF)*(nF)*(nF)*(xx)+(nF)*(nF)*(l)+(nF)*(m)+(j)]+RMi[(nF)*(nF)*(nF)*(j)+(nF)*(nF)*(l)+(nF)*(m)+(xx)])*f[nF+m]*f[nF+l];

			}
			sum3= sum3 + (1./Hi)*FMi[(2*nF)*l +i]*(-dVi[l])*FMi[(2*nF)*(m+nF) +j+nF]*f[nF+m];
			}
			sum2 = sum2+ FMi[(2*nF)*l +i]*(1./2.*dVVi[l + nF*j] + 1./2.*dVVi[j + nF*l]) ;

			sum4 = sum4 + (-3.+ep)*f[nF+i]*FMi[(2*nF)*(l+nF) +j+nF]*f[nF+l];
			
			
			}
			
            u2out[i+ j*2*nF]=0.;
            u2out[i+(j+nF)*2*nF]=0.;
			u2out[i+nF+(j)*2*nF] = (sum4 -sum2 + sum3 + 1./Hi*f[i+nF]*(-dVi[j]) +sum1 )/Hi*s ;
			u2out[i+nF+(j+nF)*2*nF]= 0.;
			
            if(i==j){
                u2out[i+nF+(j)*2*nF]=u2out[i+nF+(j)*2*nF]-1.0*(k1*k1)/(a*a)/Hi  * s ;// *a;
                u2out[i+(j+nF)*2*nF]=u2out[i+(j+nF)*2*nF] + 1./Hi  /s;// /a ;
                u2out[i+nF+(j+nF)*2*nF]= u2out[i+nF+(j+nF)*2*nF] - 3.0*Hi/Hi + ds/s/Hi ; // - 2.0*Hi/Hi;
            }

			sum1=0.0;
			sum2=0.0;
			sum3=0.0;
			sum4=0.0;

        }}

		return u2out;
	}
    
    
        // w tensor
	//
	//
	vector<double> w(vector<double> f,vector<double> p, double k1, double N)
	{
		int nT = 1;
		double a = exp(N);
		vector<double> w2out(2*nT*2*nT);
		double Hi=H(f,p);
        double s=scale(f,p,N);
        double ds=dscale(f,p,N);
		vector<double> dVVi;
		dVVi = pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);
		
        w2out[0+ 0*2*nT]=0.;
        w2out[0+1*2*nT]=+ 1./Hi;///s;
		w2out[1+(0)*2*nT] = -1.0*(k1*k1)/(a*a)/Hi;//*s ;
		w2out[1+(1)*2*nT]= - 3.0*Hi/Hi ;//+ ds/s/Hi;

		return w2out;
	}
    
		//     o                          ooooooooooo                                                           
		//    888                         88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		//   8  88         ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		//  8oooo88                           888    888          888   888          888 888     888 888        
		//o88o  o888o                        o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       

    // calculates A (the field field field term of action) for the non-trival field space metric case -- this is needed for the u2 tensor
    // A is calculated here with indices AS^{I}_{JK}

    vector<double> Acalcudd(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
		double a = exp(N);
        double Vi=pot.V(f,p);
		double Hi=H(f,p);

        vector<double> dVVi;
		dVVi=pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);
		vector<double> dVVVi;
		dVVVi=pot.dVVV(f,p);
        vector<double> Xiup(nF); 
		vector<double> A(nF*nF*nF);
        vector<double> Xid(nF);
		vector<double> fd(nF);
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		vector<double> RMi;
		RMi = fmet.Riemn(f,p); 
		vector<double> RMCi;
		RMCi = fmet.Riemncd(f,p);
      		
        double sum1=0.0;
		double sum2=0.0;
		double sum3=0.0;
		double sum4=0.0;
		double sum5=0.0;
		double sum6=0.0;
		double sumxx1=0.0;
		double s1=0.0;
		double s2=0.0;
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s2 = s2 + FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF + j];} fd[i]=s2; s2=0.0;}
		
		for(int i=0;i<nF;i++){sumxx1=sumxx1+f[nF+i]*fd[i];}

		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){
			sum1=sum1  + FMi[(2*nF)*(i)+(j)]*(-dVi[j]);
		}
			Xiup[i] = f[nF+i]/Hi*sumxx1 + 2.*(sum1 - 3.*Hi*f[nF+i]);
			sum1 =0.0;
		}//Taken into account the covariant time derivative of phi dot
			
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s1 = s1 + Xiup[j]*FMi[(2*nF)*(i+nF)+j+nF];} Xid[i]=s1; s1=0.0;}

		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			for(int l=0;l<nF;l++){
				
				sum2 =sum2 + 1./6.*FMi[(2*nF)*i+l]*(dVVVi[l + j*nF +k* nF*nF] +dVVVi[l + k*nF +j* nF*nF] + dVVVi[j + l*nF +k* nF*nF]
							+ dVVVi[j + k*nF +l* nF*nF] + dVVVi[k + l*nF +j* nF*nF] + dVVVi[k + j*nF +l* nF*nF]);				
				sum3 =sum3 + FMi[(2*nF)*i+l]* 1./2.*(dVVi[l + k*nF] + dVVi[k + l*nF]);
				sum4 =sum4 + FMi[(2*nF)*i+(l)]*1./2.*(dVVi[l + j*nF] + dVVi[j + l*nF]);//sum4 is redundant when sum3 is determined, You could optimise this by forming an array instead
				for(int m=0; m<nF; m++){for(int n=0; n<nF; n++){				
					sum5 = sum5 + 1./6. * (RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(n)+(nF)*(j)+(m)] + RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(j)+(nF)*(n)+(m)] )*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m]*fd[k]
								+ 1./6. * (RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(k)+(nF)*(n)+(m)] + RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(n)+(nF)*(k)+(m)] )*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m]*fd[j]
								+ 1./6. * (RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(j)+(nF)*(k)+(m)] + RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(k)+(nF)*(j)+(m)] )*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m]*fd[i];
				
					sum6 = sum6 + 1./6. * (  RMCi[(nF)*(nF)*(nF)*(nF)*n+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*j+k] 
										    + RMCi[(nF)*(nF)*(nF)*(nF)*k+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*j+n]
											+ RMCi[(nF)*(nF)*(nF)*(nF)*j+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*n+k] 
											+ RMCi[(nF)*(nF)*(nF)*(nF)*n+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*k+j]
											+ RMCi[(nF)*(nF)*(nF)*(nF)*k+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*n+j] 
											+ RMCi[(nF)*(nF)*(nF)*(nF)*j+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*k+n])*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m];
				
				}}
				}		
			A[i + j*nF +k* nF*nF] = -1./3. * sum2//moved to loop l
			- 1./3.*f[nF + i]/2./Hi*1./2.*( dVVi[j + k*nF] + dVVi[k + j*nF])
            - 1./3.*fd[j]/2./Hi* sum3//move to loop l 
            - 1./3.*fd[k]/2./Hi* sum4//move to loop l 
			+ 1./3.*f[nF + i] *fd[j]/8./Hi/Hi * Xid[k]
            + 1./3.*f[nF + i] * fd[k]/8./Hi/Hi * Xid[j]
            + 1./3.*fd[j] * fd[k]/8./Hi/Hi * Xiup[i]
			+ 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xid[j] *Xid[k]
            + 1./3.*fd[j]/32./Hi/Hi/Hi * Xiup[i] *Xid[k]
            + 1./3.*fd[k]/32./Hi/Hi/Hi * Xiup[i] *Xid[j]
			+ 1.*f[nF + i]*fd[j]*fd[k]/8./Hi/Hi/Hi*2.*Vi
			- 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xid[j] * Xid[k] * (k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            - 1./3.*fd[j]/32./Hi/Hi/Hi * Xiup[i] * Xid[k] * (k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.
            - 1./3.*fd[k]/32./Hi/Hi/Hi * Xiup[i] * Xid[j] * (k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4.;
    		A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] + 1./3.*f[nF+i]*FMi[(2*nF)*(j+nF)+k+nF]/2./Hi*(-k2*k2-k3*k3+k1*k1)/a/a/2.;
			if(i==k){A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] + 1./3.*fd[j]/2./Hi*(-k1*k1-k3*k3+k2*k2)/a/a/2.;}
			if(i==j){A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] + 1./3.*fd[k]/2./Hi*(-k2*k2-k1*k1+k3*k3)/a/a/2.;}
			
			A[i + j*nF +k* nF*nF] = A[i + j*nF +k* nF*nF] - 1./2.* sum5/Hi + 1./3. * sum6;

			sum2=0.0;
			sum3=0.0;
			sum4=0.0;
			sum5=0.0;
			sum6=0.0;
            }}}

        
        return A;
    }	

	
		//     o       oooooooo8                ooooooooooo                                                           
		//    888     888                       88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		//   8  88     888oooooo ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		//  8oooo88           888                   888    888          888   888          888 888     888 888        
		//o88o  o888o o88oooo888                   o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       

    //calculates AS (the "slow" parts of the field field field term of action for non-trival field space case
    //this is used only for initial conditions and is calculated with indices AS^{I}_{JK}
    vector<double> AScalcudd(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
        double Vi=pot.V(f,p);
		double Hi=H(f,p);
     
        vector<double> dVVi;
		dVVi=pot.dVV(f,p);
		vector<double> dVi;
		dVi =  pot.dV(f,p);
		vector<double> dVVVi;
		dVVVi=pot.dVVV(f,p);
        vector<double> Xiup(nF); 
		vector<double> AS(nF*nF*nF);
        vector<double> Xid(nF);
		vector<double> fd(nF);
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		vector<double> RMi;
		RMi = fmet.Riemn(f,p); 
		vector<double> RMCi;
		RMCi = fmet.Riemncd(f,p);

        double sum1=0.0;
		double sum2=0.0;
		double sum3=0.0;
		double sum4=0.0;
		double sum5=0.0;
		double sum6=0.0;
		double sumxx1=0.0;
		
		double s1=0.0;
		double s2=0.0;
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s2 = s2 + FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF + j];}fd[i]=s2; s2=0.0;}
		
		//for(int i=0;i<nF;i++){sum1=sum1+fd[i]*f[nF+i];}
		
		for(int i=0;i<nF;i++){sumxx1=sumxx1+f[nF+i]*fd[i];}

		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){
			sum1=sum1  + FMi[(2*nF)*(i)+(j)]*(-dVi[j]);
		}
			Xiup[i] = f[nF+i]/Hi*sumxx1 + 2.*(sum1 - 3.*Hi*f[nF+i]);
			sum1 =0.0;
		}//Taken into account the covariant time derivative of phi dot
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s1 = s1 + Xiup[j]*FMi[(2*nF)*(i+nF)+j+nF];} Xid[i]=s1; s1=0.0;}

		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			for(int l=0;l<nF;l++){
				
				sum2 =sum2 + 1./6.*FMi[(2*nF)*i+l]*(dVVVi[l + j*nF +k* nF*nF] +dVVVi[l + k*nF +j* nF*nF] + dVVVi[j + l*nF +k* nF*nF]
							+ dVVVi[j + k*nF +l* nF*nF] + dVVVi[k + l*nF +j* nF*nF] + dVVVi[k + j*nF +l* nF*nF]);
				sum3 =sum3 + FMi[(2*nF)*i+l]* 1./2.*(dVVi[l + k*nF] + dVVi[k + l*nF]);
				sum4 =sum4 + FMi[(2*nF)*i+(l)]*1./2.*(dVVi[l + j*nF] + dVVi[j + l*nF]);///sum4 is redundant when sum3 is determined, You could optimise this by forming an array instead
				for(int m=0; m<nF; m++){for(int n=0; n<nF; n++){				
					sum5 = sum5 + 1./6. * (RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(n)+(nF)*(j)+(m)] + RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(j)+(nF)*(n)+(m)] )*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m]*fd[k]
								+ 1./6. * (RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(k)+(nF)*(n)+(m)] + RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(n)+(nF)*(k)+(m)] )*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m]*fd[j]
								+ 1./6. * (RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(j)+(nF)*(k)+(m)] + RMi[(nF)*(nF)*(nF)*(l)+(nF)*(nF)*(k)+(nF)*(j)+(m)] )*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m]*fd[i];
								
					sum6 = sum6 + 1./6. * (  RMCi[(nF)*(nF)*(nF)*(nF)*n+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*j+k] 
										    + RMCi[(nF)*(nF)*(nF)*(nF)*k+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*j+n]
											+ RMCi[(nF)*(nF)*(nF)*(nF)*j+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*n+k] 
											+ RMCi[(nF)*(nF)*(nF)*(nF)*n+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*k+j]
											+ RMCi[(nF)*(nF)*(nF)*(nF)*k+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*n+j] 
											+ RMCi[(nF)*(nF)*(nF)*(nF)*j+(nF)*(nF)*(nF)*l+(nF)*(nF)*m+(nF)*k+n])*FMi[(2*nF)*(i)+n]*f[nF+l]*f[nF+m];
				
				}}
				}		
			AS[i + j*nF +k* nF*nF] = -1./3. * sum2//moved to loop l
			- 1./3.*f[nF + i]/2./Hi*1./2.*( dVVi[j + k*nF] + dVVi[k + j*nF])
            - 1./3.*fd[j]/2./Hi* sum3//move to loop l 
            - 1./3.*fd[k]/2./Hi* sum4//move to loop l 
			+ 1./3.*f[nF + i] *fd[j]/8./Hi/Hi * Xid[k]
            + 1./3.*f[nF + i] * fd[k]/8./Hi/Hi * Xid[j]
            + 1./3.*fd[j] * fd[k]/8./Hi/Hi * Xiup[i]
			+ 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xid[j] *Xid[k]
            + 1./3.*fd[j]/32./Hi/Hi/Hi * Xiup[i] *Xid[k]
            + 1./3.*fd[k]/32./Hi/Hi/Hi * Xiup[i] *Xid[j]
			+ 1.*f[nF + i]*fd[j]*fd[k]/8./Hi/Hi/Hi*2.*Vi
			- 1./3.*f[nF + i]/32./Hi/Hi/Hi * Xid[j] * Xid[k] * (k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            - 1./3.*fd[j]/32./Hi/Hi/Hi * Xiup[i] * Xid[k] * (k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.
            - 1./3.*fd[k]/32./Hi/Hi/Hi * Xiup[i] * Xid[j] * (k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4.;
			
			AS[i + j*nF +k* nF*nF] = AS[i + j*nF +k* nF*nF] - 1./2.* sum5/Hi + 1./3. * sum6;
			sum2=0.0;
			sum3=0.0;
			sum4=0.0;
			sum5=0.0;
			sum6=0.0;
            }}}

        return AS;
    }	
	
		//oooooooooo                      ooooooooooo                                                           
		// 888    888                     88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		// 888oooo88       ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		// 888    888                         888    888          888   888          888 888     888 888        
		//o888ooo888                         o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       

	//Calculates B term of action for non-trivial field space metric with indices B^{I}_{JK}
    vector<double> Bcalcudd(vector<double> f,vector<double> p, double k1, double k2, double k3,double N)
	{
		
        double Hi=H(f,p);
        
		vector<double> dVVi;
       // dVVi = new double[nF*nF];
		dVVi=pot.dVV(f,p);
		vector<double> dVi; //dVi = new double[nF];
		dVi =  pot.dV(f,p);
		vector<double> dVVVi; //dVVVi = new double[nF*nF*nF];
		dVVVi=pot.dVVV(f,p);
        vector<double> Xiup(nF);
		vector<double> Xid(nF);
		vector<double> B(nF*nF*nF);
		vector<double> fd(nF);
        vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		vector<double> RMi;
		RMi = fmet.Riemn(f,p); 
		
        double sum1=0;
		double sum2=0;
		double s1=0.0;
		double s2=0.0;
		double sumxx1=0.0;
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s2 = s2 + FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF + j];}fd[i]=s2; s2=0.0;}
		
		//for(int i=0;i<nF;i++){sum1=sum1+fd[i]*f[nF+i];}			
		
		for(int i=0;i<nF;i++){sumxx1=sumxx1+f[nF+i]*fd[i];}

		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){
			sum1=sum1  + FMi[(2*nF)*(i)+(j)]*(-dVi[j]);
		}
			Xiup[i] = f[nF+i]/Hi*sumxx1 + 2.*(sum1 - 3.*Hi*f[nF+i]);
			sum1 =0.0;
		}//Taken into account the covariant time derivative of phi dot
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s1 = s1 + Xiup[j]*FMi[(2*nF)*(i+nF)+j+nF];} Xid[i]=s1; s1=0.0;}
		
        
        for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			for(int l=0;l<nF;l++){for(int n=0;n<nF;n++){
				sum2 = sum2 + 1./2. *( RMi[(nF)*(nF)*(nF)*(k)+(nF)*(nF)*(n)+(nF)*(j)+(l)] +RMi[(nF)*(nF)*(nF)*(k)+(nF)*(nF)*(j)+(nF)*(n)+(l)] 
			)*FMi[(2*nF)*(n)+(i)]*f[nF+l];// David fixed the index notation here, symmetrised over ij only as the kth index is the momentum so only 2 Riemann tensors and hence factor of 1/2
			}}
			B[i + j*nF +k* nF*nF] = 1.*f[nF + i]*fd[j]*fd[k]/4./Hi/Hi
			- 1./2.*f[nF + i] * fd[k]/8./Hi/Hi/Hi * Xid[j]
            - 1./2.*fd[j] * fd[k]/8./Hi/Hi/Hi * Xiup[i]
			+ 1./2.*f[nF + i] * fd[k]/8./Hi/Hi/Hi * Xid[j]*(k2*k2+k3*k3 - k1*k1)*(k2*k2+k3*k3 - k1*k1)/k2/k2/k3/k3/4.
            + 1./2.*fd[j] * fd[k]/8./Hi/Hi/Hi * Xiup[i]*(k1*k1+k3*k3 - k2*k2)*(k1*k1+k3*k3 - k2*k2)/k1/k1/k3/k3/4.;
			B[i + j*nF +k* nF*nF] = B[i + j*nF +k* nF*nF] - FMi[(2*nF)*(j+nF)+(k+nF)]*1.*Xiup[i]/4./Hi*(-k1*k1-k2*k2+k3*k3)/k1/k1/2.;
			B[i + j*nF +k* nF*nF] = B[i + j*nF +k* nF*nF] + 4./3.*sum2;
			if(i==k){B[i + j*nF +k* nF*nF] = B[i + j*nF +k* nF*nF] - 1.*Xid[j]/4./Hi*(-k1*k1-k2*k2+k3*k3)/k2/k2/2.;}
			sum2 = 0.0;
			
			
		}}}
		
			
        return B;
    }

		//  oooooooo8                     ooooooooooo                                                           
		//o888     88                     88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		//888              ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		//888o     oo                         888    888          888   888          888 888     888 888        
		// 888oooo88                         o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       
	
	
	//Calculates C term of action for non-trivial field space metric case with indices C^{I}_{JK}
    vector<double> Ccalcudd(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{	
		double Hi=H(f,p);
        vector<double> FMi;
		FMi = fmet.fmetric(f,p);
     	vector<double> dVVi; //dVVi = new double[nF*nF];
		dVVi=pot.dVV(f,p);
		vector<double> dVi; //dVi = new double[nF];
		dVi =  pot.dV(f,p);
		vector<double> dVVVi; //dVVVi = new double[nF*nF*nF];
		dVVVi=pot.dVVV(f,p);
        vector<double> Xid(nF); vector<double> C(nF*nF*nF);
		vector<double> Xiup(nF);
		vector<double> fd(nF);
		vector<double> CHR;
		CHR = fmet.Chroff(f,p);
		vector<double> u1 =u(f,p);
        double sum1=0;
		double s1=0.0;
		double s2=0.0;
		double sumxx1=0.0;
	
		for(int i=0;i<nF;i++){
			for(int j=0;j<nF;j++){
			s2 = s2 + FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF + j];}
			fd[i]=s2; 
			s2=0.0;}
		//for(int i=0;i<nF;i++){sum1=sum1+fd[i]*f[nF+i];}	
		for(int i=0;i<nF;i++){sumxx1=sumxx1+f[nF+i]*fd[i];}

		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){
			sum1=sum1  + FMi[(2*nF)*(i)+(j)]*(-dVi[j]);
		}
			Xiup[i] = f[nF+i]/Hi*sumxx1 + 2.*(sum1 - 3.*Hi*f[nF+i]);
			sum1 =0.0;
			
		}//Taken into account the covariant time derivative of phi dot
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s1 = s1 + Xiup[j]*FMi[(2*nF)*(i+nF)+j+nF];} Xid[i]=s1; s1=0.0;}
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			C[i + j*nF +k* nF*nF] = 1.*f[nF + i]*fd[j]*fd[k]/8./Hi/Hi/Hi
			- 1.*f[nF + i] * fd[j] *fd[k]/8./Hi/Hi/Hi *(k1*k1+k2*k2 - k3*k3)*(k1*k1+k2*k2 - k3*k3)/k1/k1/k2/k2/4. ;
			if(i==j){C[i + j*nF +k* nF*nF] = C[i + j*nF +k* nF*nF] - 1.*fd[k]/2./Hi;}
			C[i + j*nF +k* nF*nF] = C[i + j*nF +k* nF*nF] + FMi[(2*nF)*(j+nF)+(k+nF)]*f[nF+i]/2./Hi*(-k1*k1-k3*k3+k2*k2)/k1/k1/2.;
			if(i==k){C[i + j*nF +k* nF*nF] = C[i + j*nF +k* nF*nF] + fd[j]/2./Hi*(-k2*k2-k3*k3+k1*k1)/k2/k2/2.;}

		}}}
		
        return C;
    }
	
    
    
    
		//     o      oooooooooo    oooooooo8                     oooooooooo             ooooo                   oooo                        
		//    888      888    888 o888     88                      888    888 ooooooooo8  888  oo oooooo    ooooo888  ooooooooo8 oooo   oooo 
		//   8  88     888oooo88  888              ooooooooo       888oooo88 888oooooo8   888   888   888 888    888 888oooooo8    888o888   
		//  8oooo88    888    888 888o     oo                      888  88o  888          888   888   888 888    888 888           o88 88o   
		//o88o  o888o o888ooo888   888oooo88                      o888o  88o8  88oooo888 o888o o888o o888o  88ooo888o  88oooo888 o88o   o88o 
	
    
    // Rearanges the indices of the A tensor to A^{IJK}
	vector<double> Acalcuuu(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	 {
		vector<double> A;

		vector<double> Aout(nF*nF*nF);
		
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		A = Acalcudd(f,p, k1, k2, k3, N);
		double sum =0.0;
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){for(int l=0;l<nF;l++){for(int m=0;m<nF;m++){
			sum = sum + FMi[(2*nF)*(m)+k]*FMi[(2*nF)*(l)+j]*A[i+nF*l+nF*nF*m];
		}}
			Aout[i+nF*j+nF*nF*k]=sum;
	
			sum=0.0;
		}}}
		return Aout;
	 }
	 
	 // Rearanges the indices of the AS tensor to AS^{IJK}
	 vector<double> AScalcuuu(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	 {
		vector<double> AS;

		vector<double> ASout(nF*nF*nF);
		
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		AS = AScalcudd(f,p, k1, k2, k3, N);
		double sum =0.0;
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){for(int l=0;l<nF;l++){for(int m=0;m<nF;m++){
			sum = sum + FMi[(2*nF)*(m)+k]*FMi[(2*nF)*(l)+j]*AS[i+nF*l+nF*nF*m];
		}}
			ASout[i+nF*j+nF*nF*k]=sum;

			sum=0.0;
		}}}
		return ASout;
	 }
	 
	 // Rearanges the indices of the B tensor to B^{IJK}
	 vector<double> Bcalcuuu(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	 {
		vector<double> B;

		vector<double> Bout(nF*nF*nF);
		
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		B = Bcalcudd(f,p, k1, k2, k3, N);
		double sum =0.0;
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){for(int l=0;l<nF;l++){for(int m=0;m<nF;m++){
			sum = sum + FMi[(2*nF)*(m)+k]*FMi[(2*nF)*(l)+j]*B[i+nF*l+nF*nF*m];
		}}
			Bout[i+nF*j+nF*nF*k]=sum;

			sum=0.0;
		}}}
		return Bout;
	 }
	 
	 // Rearanges the indices of the A tensor to B_{IJ}^{K}
	 vector<double> Bcalcddu(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	 {
		vector<double> B;

		vector<double> Bout(nF*nF*nF);
		
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		B = Bcalcudd(f,p, k1, k2, k3, N);
		double sum =0.0;
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){for(int l=0;l<nF;l++){for(int m=0;m<nF;m++){
			sum = sum + FMi[(2*nF)*(m+nF)+i+nF]*FMi[(2*nF)*(k)+l]*B[m+nF*j+nF*nF*l];
		}}
			Bout[i+nF*j+nF*nF*k]=sum;

			sum=0.0;
		}}}
		return Bout;
	 }
	
	// Rearanges the indices of the C tensor to C^{IJK}
	vector<double> Ccalcuuu(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	 {
		vector<double> C;

		vector<double> Cout(nF*nF*nF);
		
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		C = Ccalcudd(f,p, k1, k2, k3, N);
		double sum =0.0;
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){for(int l=0;l<nF;l++){for(int m=0;m<nF;m++){
			sum = sum + FMi[(2*nF)*(m)+k]*FMi[(2*nF)*(l)+j]*C[i+nF*l+nF*nF*m];
		}}
			Cout[i+nF*j+nF*nF*k]=sum;

			
			sum=0.0;
		}}}
		return Cout;
	 }

	 
	 // Rearanges the indices of the C tensor to C_{IJ}^K
	 vector<double> Ccalcddu(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	 {
		vector<double> C;

		vector<double> Cout(nF*nF*nF);
		
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		
		C = Ccalcudd(f,p, k1, k2, k3, N);
		double sum =0.0;
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){for(int l=0;l<nF;l++){for(int m=0;m<nF;m++){
			sum = sum + FMi[(2*nF)*(i+nF)+m+nF]*FMi[(2*nF)*(k)+l]*C[m+nF*j+nF*nF*l];
		}}
			Cout[i+nF*j+nF*nF*k]=sum;

			sum=0.0;
		}}}
		return Cout;
	 }
	 
	 
		//ooooo  oooo  ooooooo                       ooooooooooo                                                           
		// 888    88 o88    888o                     88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		// 888    88     88888o       ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		// 888    88 88o    o888                         888    888          888   888          888 888     888 888        
		//  888oo88    88ooo88                          o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       

    
	//calculates u3
	vector<double> u(vector<double> f, vector<double> p, double k1, double k2, double k3,double N)
	{
        vector<double>  A, B,B2, B3, C, C2,C3;
        double Hi;
		Hi=H(f,p);
        double s=scale(f,p,N);
        
        A = Acalcudd(f,p, k1, k2, k3 ,N);
        B= Bcalcddu(f,p, k2, k3, k1 ,N);
        B2=  Bcalcudd(f,p, k1, k2, k3 ,N);
        B3=Bcalcudd(f,p, k1, k3, k2 ,N);
        C=  Ccalcudd(f,p, k1, k2, k3 ,N);
        C2=  Ccalcudd(f,p, k1, k3, k2 ,N);
        C3 = Ccalcddu(f,p, k3, k2, k1 ,N);


        vector<double> u3out(2*nF*2*nF*2*nF);
		
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){for(int k=0;k<nF;k++){
			u3out[i+j*2*nF+k*2*nF*2*nF]= -B[j+k*nF+i*nF*nF]/Hi;
			
            u3out[(i)+(nF+j)*2*nF+k*2*nF*2*nF]= -C[i+j*nF+k*nF*nF]/Hi  /s; // /a;
            u3out[(i)+j*2*nF+(k+nF)*2*nF*2*nF]= -C2[i+k*nF+j*nF*nF]/Hi /s;// /a;
			
			u3out[(i)+(j+nF)*2*nF+(k+nF)*2*nF*2*nF]= 0.;
			
            u3out[(nF+i) + j*2*nF + k*2*nF*2*nF]= 3.*A[i+j*nF+k*nF*nF]/Hi  *s;// *a;
		
			u3out[(nF+i)+(nF+j)*2*nF+k*2*nF*2*nF]=B3[i+k*nF+j*nF*nF]/Hi ;//momentum
			u3out[(nF+i)+(j)*2*nF+(k+nF)*2*nF*2*nF]=B2[i+j*nF+k*nF*nF]/Hi ;//field
			
            u3out[(nF+i)+(j+nF)*2*nF + (k+nF)*2*nF*2*nF]=C3[k+j*nF+i*nF*nF]/Hi  /s;// /a;k down i up


		}}}
        return u3out;
	}

		//oooo   oooo  oo                      ooooooooooo                                                           
		// 8888o  88 o888                      88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		// 88 888o88  888       ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		// 88   8888  888                          888    888          888   888          888 888     888 888        
		//o88o    88 o888o                        o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       

	//calculates N1_I tesnor
	vector<double> N1(vector<double> f,vector<double> p, double N)
	{
		double Hi=H(f,p);
		vector<double> dVi;
		vector<double> Ni(2*nF);
		double ep=Ep(f,p);
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		dVi=pot.dV(f,p);
		double sum=0.0;
		for(int i=0;i<nF;i++){
			sum=0.0;
			for(int l=0;l<nF;l++){
			sum =  sum +  FMi[(2*nF)*(l+nF)+i+nF]*f[nF+l];
			}
			Ni[i]=-1./2.0/Hi/ep*sum;
			Ni[nF+i] = 0. ;
			
			}
			
		return Ni;
	}

		//oooo   oooo  ooooooo                       ooooooooooo                                                           
		// 8888o  88 o88     888                     88  888  88 ooooooooo8 oo oooooo    oooooooo8    ooooooo  oo oooooo   
		// 88 888o88       o888       ooooooooo          888    888oooooo8   888   888  888ooooooo  888     888 888    888 
		// 88   8888    o888   o                         888    888          888   888          888 888     888 888        
		//o88o    88 o8888oooo88                        o888o     88oooo888 o888o o888o 88oooooo88    88ooo88  o888o       


	 // Calculates the N2_{IJ} tensor
	vector<double> N2(vector<double> f, vector<double> p, double k1, double k2, double k3, double N)
	{
		double Hd=Hdot(f,p);
		double Hin=H(f,p);
		vector<double> dVi, dVVi;
		vector<double> Nii(2*nF*2*nF);
		double s = scale(f,p,N);
		dVi=pot.dV(f,p);
		dVVi=pot.dVV(f,p);
		vector<double> fd(nF);
		vector<double> N11;
		N11= N1(f,p,N);
		vector<double> FMi;
		FMi = fmet.fmetric(f,p);
		vector<double> RMi;
		RMi = fmet.Riemn(f,p); 
		vector<double> CHi;
		CHi = fmet.Chroff(f,p);

		double sum3 = 0.0;


		for(int i=0;i<nF;i++){sum3=sum3+dVi[i]*f[nF+i]/Hin/Hin/Hin;}
		
		double s2=0.0;
		for(int i=0;i<nF;i++){for(int j=0;j<nF;j++){s2 = s2 + FMi[(2*nF)*(i+nF)+(j+nF)]*f[nF + j];} fd[i]=s2; s2=0.0;}

		double ep = -Hd/Hin/Hin;
		for(int i=0;i<nF;i++){for(int j=0; j<nF; j++){
		Nii[i+(j)*2*nF]= 2./ep/Hin/Hin/6. * (fd[i]*fd[j] *(-3./2. + 9./2./ep + 3./4.*sum3/ep/ep));
		Nii[i+(j+nF)*2*nF]=2./ep/Hin/Hin/6.*3./2.*fd[i]*fd[j]/Hin/ep /s;
		Nii[i+nF+(j)*2*nF]=2./ep/Hin/Hin/6.*3./2.*fd[i]*fd[j]/Hin/ep /s;
		Nii[i+nF+(j+nF)*2*nF]=0.;
		Nii[i+nF+(j)*2*nF] = Nii[i+nF + (j) * 2*nF] - FMi[(2*nF)*(i+nF)+j+nF]*2./ep/Hin/Hin/6. * 3./2.*Hin/k1/k1*((-k2*k2-k3*k3+k1*k1)/2. + k3*k3)  /s;// /a ;
		Nii[i+(j+nF)*2*nF] = Nii[i + (j+nF) * 2*nF] - FMi[(2*nF)*(i+nF)+j+nF]*2./ep/Hin/Hin/6. * 3./2.*Hin/k1/k1*((-k2*k2-k3*k3+k1*k1)/2. + k2*k2)  /s;// /a;}
		
		

		}}
		
		
		
		return Nii;
		
	}

};
#endif
