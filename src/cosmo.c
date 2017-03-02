//==========================================================================//
//  Copyright (c) 2013       Cullan Howlett & Marc Manera,                  //
//                           Institute of Cosmology and Gravitation,        //
//                           University of Portsmouth.                      //
//                                                                          //
//  MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017           //
//                                                                          //
//  This file is part of PICOLA.                                            //
//                                                                          //
//  PICOLA is free software: you can redistribute it and/or modify          //
//  it under the terms of the GNU General Public License as published by    //
//  the Free Software Foundation, either version 3 of the License, or       //
//  (at your option) any later version.                                     //
//                                                                          //
//  PICOLA is distributed in the hope that it will be useful,               //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//  GNU General Public License for more details.                            //
//                                                                          //
//  You should have received a copy of the GNU General Public License       //
//  along with PICOLA.  If not, see <http://www.gnu.org/licenses/>.         //
//==========================================================================//

//==========================================================================//
// This file contains all the cosmology related routines                    //
//==========================================================================//

#include "vars.h"
#include "proto.h"

//=============================================================
// This file contains updated cosmo function 
// The growth-factor routines found at the end of this file are 
// not (should not) be used in the main code
//=============================================================
#include "new_cosmo.h"

//=============================================================
// Functions for COLA modified time-stepping (used when StdDA=0)
// This is updated to work with a general cosmology
//=============================================================
double Sq(double ai,double af,double aRef) {
  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

  F.function = &fun;
  F.params = &alpha;

  gsl_integration_qag(&F, (double)ai, (double)af, 0, 1e-5, 5000, 6, w, &result, &error); 
  gsl_integration_workspace_free (w);

  if (fullT == 1) {
    return result/gpQ(aRef);
  } else {
    return result;
  }
}

//=============================================================
// The integrand for integration above
//=============================================================
double fun(double a, void * params) {   
  double f;
  if (fullT == 1) {
    f = gpQ(a)/Qfactor(a); 
  } else {
    f = 1.0/Qfactor(a);
  }   
  return f;
}

double Sphi(double ai,double af,double aRef) {
  return (gpQ(af)-gpQ(ai))*aRef/Qfactor(aRef)/DERgpQ(aRef);
}

double gpQ(double a) { 
  return pow(a,nLPT);      
}

//=============================================================
// This must return d(gpQ)/da
//=============================================================
double DERgpQ(double a) { 
  return nLPT * pow(a,nLPT-1.0);
}

//============================================================
// Functions for Quinn et al time-stepping (used when StdDA=2)
//============================================================
double SqStd(double ai,double af) {
  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

  F.function = &funSqStd;
  F.params = &alpha;

  gsl_integration_qag(&F,ai,af,0,1e-5,5000,6,w,&result,&error); 
  gsl_integration_workspace_free (w);

  return result;
}

double funSqStd(double a, void * params) {
  return 1.0/Qfactor(a);
}

double SphiStd(double ai,double af) {
  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

  F.function = &funSphiStd;
  F.params = &alpha;

  gsl_integration_qag (&F,ai,af,0,1e-5,5000,6,w,&result,&error); 
  gsl_integration_workspace_free (w);
     
  return result;
}

double funSphiStd(double a, void * params) {       
  return a/Qfactor(a);
}
 
//=================================================
// The functions below are used when stepDistr == 2
//=================================================

// Holds parameters needed for the integration
struct quadratic_params { 
  double y;
};

//===============================
// Solves y = CosmoTime(a) for a
//===============================
double AofTime(double y) {
  int status;
  int iter = 0, max_iter = 1000;
  double r = 0.0;
  double x_lo = 0.0001, x_hi = 2.01;
  struct quadratic_params params = {y};

  gsl_function F;
  gsl_root_fsolver *s;
  const gsl_root_fsolver_type * T;

  F.function = &AofTimeFun;
  F.params = &params;
     
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
     
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo,x_hi,0,1.0e-6);
  } while (status == GSL_CONTINUE && iter < max_iter);
     
  gsl_root_fsolver_free (s);
     
  return r;
}

double AofTimeFun(double a,void * params) {
  struct quadratic_params * p = (struct quadratic_params *) params;
  double y = p->y;
  return CosmoTime(a)-y;
}
  
//=============================================================     
// As defined below, CosmoTime(a)=int_0^a CosmoTimeFun(a) da
// which corresponds to the usual time coordinate scaled by H0.
//=============================================================     
double CosmoTime(double af) {
  double ai = 1.0e-8;
  double alpha = 0.0;
  double result, error;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
     
  F.function = &CosmoTimeFun;
  F.params = &alpha;
     
  gsl_integration_qag (&F,ai,af,0,1e-5,5000,6,w, &result, &error); 
    
  gsl_integration_workspace_free (w);
     
  return result;
}

double CosmoTimeFun (double a, void * params) { 
  return a*a/Qfactor(a);
}

//===========================================================================
//===========================================================================
// 
// Deprecated methods from the original PICOLA code
// These have been replaced. Only used in new_cosmo.h to
// compare to the new version
//
//===========================================================================
//===========================================================================

//===============================
// Q(a) = a^3 H(a)/H0 for LCDM
//===============================
double Qfactor_LCDM(double a) { 
  return sqrt(Omega/(a*a*a)+(1.0-Omega))*a*a*a;
}

//===================================
// First order Growth factor for LCDM
// normalized to unity at a = 1 
//===================================
double growthD(double a){ 
    return growthDtemp(a)/growthDtemp(1.0);
}

//===================================
// Fit to first order growth-factor
//===================================
double growthDtemp(double a) {
  // Decided to use the analytic expression for LCDM. More transparent if I change this to numerical integration? YES!
  double hyperP = 0.0, hyperM = 0.0;
  double x = -Omega/(Omega-1.0)/(a*a*a);
    
  if (fabs(x-1.0) < 1.e-3) {
    hyperP = 0.859596768064608  - 0.1016599912520404*(-1.0 + x)  + 0.025791094277821357*pow(-1.0 + x,2) 
           - 0.008194025861121475*pow(-1.0 + x,3) + 0.0029076305993447644*pow(-1.0 + x,4) 
           - 0.0011025426387159761*pow(-1.0 + x,5) + 0.00043707304964624546*pow(-1.0 + x,6) 
           - 0.0001788889964687831*pow(-1.0 + x,7);
    hyperM = 1.1765206505266006 + 0.15846194123099624*(-1.0 + x) - 0.014200487494738975*pow(-1.0 + x,2) 
           + 0.002801728034399257*pow(-1.0 + x,3) - 0.0007268267888593511*pow(-1.0 + x,4) 
           + 0.00021801569226706922*pow(-1.0 + x,5) - 0.00007163321597397065*pow(-1.0 + x,6) 
           + 0.000025063737576245116*pow(-1.0 + x,7);
  } else {
    if (x < 1.0) {
      hyperP = gsl_sf_hyperg_2F1(1.0/2.0,2.0/3.0,5.0/3.0,-x);
      hyperM = gsl_sf_hyperg_2F1(-1.0/2.0,2.0/3.0,5.0/3.0,-x);
    }
    x=1.0/x;
    if ((x < 1.0) && (x>1.0/30.0)) {
      hyperP  = gsl_sf_hyperg_2F1(-1.0/6.0,0.5,5.0/6.0,-x);
      hyperP *= 4*sqrt(x);
      hyperP += -3.4494794123063873799*pow(x,2.0/3.0);
      
      hyperM  = gsl_sf_hyperg_2F1(-7.0/6.0,-0.5,-1.0/6.0,-x);
      hyperM *= 4.0/7.0/sqrt(x);
      hyperM += pow(x,2.0/3.0)*(-1.4783483195598803057); //-(Gamma[-7/6]*Gamma[5/3])/(2*sqrt[Pi])
    }
    if (x<=1.0/30.0) {
      hyperP = 4.0*sqrt(x) - 3.4494794123063865*pow(x,2.0/3.0) + 0.4*pow(x,1.5) - 0.13636363636363635*pow(x,2.5) 
             + 0.07352941176470587*pow(x,3.5) - 0.04755434782608695*pow(x,4.5) + 0.033943965517241374*pow(x,5.5) 
             - 0.02578125*pow(x,6.5) + 0.020436356707317072*pow(x,7.5) - 0.01671324384973404*pow(x,8.5) 
             + 0.013997779702240564*pow(x,9.5) - 0.011945562847590041*pow(x,10.5) + 0.01035003662109375*pow(x,11.5) 
             - 0.009080577904069926*pow(x,12.5);
      hyperM = 0.5714285714285715/sqrt(x) + 2.0*sqrt(x) - 1.4783483195598794*pow(x,2.0/3.0) + 0.1*pow(x,1.5) 
             - 0.022727272727272735*pow(x,2.5) + 0.009191176470588237*pow(x,3.5) - 0.004755434782608697*pow(x,4.5) 
             + 0.002828663793103449*pow(x,5.5) - 0.0018415178571428578*pow(x,6.5) + 0.0012772722942073172*pow(x,7.5) 
             - 0.0009285135472074472*pow(x,8.5) + 0.0006998889851120285*pow(x,9.5) - 0.0005429801294359111*pow(x,10.5) 
             + 0.0004312515258789064*pow(x,11.5) - 0.00034925299631038194*pow(x,12.5);
    }
  }
    
  if (a>0.2) {
    return sqrt(1.0 + (-1.0 + pow(a,-3))*Omega)*(3.4494794123063873799*pow(-1.0 + 1.0/Omega,2.0/3.0) 
           + (hyperP*(4*pow(a,3)*(-1.0 + Omega) - Omega) - 7.0*pow(a,3)*hyperM*(-1.0 + Omega))/(pow(a,5)*(-1.0+ Omega) - pow(a,2)*Omega));
  } else {
    return (a*pow(1 - Omega,1.5)*(1291467969*pow(a,12)*pow(-1 + Omega,4) 
            + 1956769650*pow(a,9)*pow(-1 + Omega,3)*Omega + 8000000000*pow(a,3)*(-1 + Omega)*pow(Omega,3) 
            + 37490640625*pow(Omega,4)))/(1.5625e10*pow(Omega,5));
  } 
}

//====================================
// Second order growth factor in LCDM
//====================================
double growthD2(double a) { 
  return growthD2temp(a)/growthD2temp(1.0);
}

//===========================================================
// Fitting formula for second order growth factor in LCDM
//===========================================================
double growthD2temp(double a) {
  double d = growthD(a);
  double Om = Omega/(Omega+(1.0-Omega)*a*a*a);
  return (d*d*pow(Om,-1.0/143.0));
}
 
//=========================
// T[D_{2lpt}]=dD_{2lpt}/dy
//=========================
double growthD2v(double a) { 
    double d2 = growthD2(a);
    double Om = Omega/(Omega+(1.0-Omega)*a*a*a);
    return Qfactor_LCDM(a)*(d2/a)*2.0*pow(Om,6./11.);
}

//=========================
// D_{-}, the decaying mode
//=========================
double decayD(double a) { 
  return sqrt(Omega/(a*a*a)+1.0-Omega);
}

//=======================================================================================
// Returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor_LCDM(a). 
// Note however in this code we only interested in the case where nGrowth=1.0, nDecay=0.0 
//=======================================================================================
double DprimeQ(double a) {
  double Nn = 6.0*pow(1.0 - Omega,1.5)/growthDtemp(1.0);
  return (1.0/decayD(a))*(Nn-(3.0*Omega*growthD(a))/(2.0*a));      
}

