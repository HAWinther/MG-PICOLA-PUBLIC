//==========================================================================//
//                                                                          //
//  MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017           //
//                                                                          //
// This file contains updated cosmology functions needed to solve the       //
// general growth ODEs plus look-up functions                               //
//                                                                          //
// The original functions used in PICOLA are found as                       //
// [ xxxx_LCDMFit(double a) ]                                               //
// The growth-factor routines found at the end of this file are             //
// not (should not) be used in the main code                                //
//                                                                          //
//  Below d/dy = Q(a) d/da                                                   //
//                                                                          //
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

#include "vars.h"
#include "proto.h"
#include "user_defined_functions.h"

// Accuracy parameters for numerical integration
#define MY_GSL_HSTART  1.0e-7
#define MY_GSL_EPS     1.0e-7
#define MY_GSL_REL     0.0

// Accuracy parameters for lower accuracy (makes it faster, used for second order kernel)
#define MY_GSL_HSTART_LOWACC  1.0e-4
#define MY_GSL_EPS_LOWACC     1.0e-4
#define MY_GSL_REL_LOWACC     0.0

// The scale-factor we normalize the growth-factor at
// Side-effects of changing this has not been tested
#define GROWTH_NORMALIZATION_SCALE_FACTOR 1.0

// Redshift for which to start and end the integration
#define REDSHIFT_START_INTEGRATION (MAX(200.0, Init_Redshift))
#define REDSHIFT_END_INTEGRATION   (-0.5)

//========================================================
// Splines needed to store growth-factors and derivatives
//========================================================

struct MySplineContainer{

#ifdef SCALEDEPENDENT
  GSL_2D_Spline D_of_scale_spline;
  GSL_2D_Spline dDdy_of_scale_spline;
  GSL_2D_Spline ddDddy_of_scale_spline;
  GSL_2D_Spline D2_of_scale_spline;
  GSL_2D_Spline dD2dy_of_scale_spline;
  GSL_2D_Spline ddD2ddy_of_scale_spline;
#ifdef MASSIVE_NEUTRINOS
  GSL_2D_Spline DLCDM_of_scale_spline;
  GSL_2D_Spline dDLCDMdy_of_scale_spline;
  GSL_2D_Spline ddDLCDMddy_of_scale_spline;
  GSL_2D_Spline D2LCDM_of_scale_spline;
  GSL_2D_Spline dD2LCDMdy_of_scale_spline;
  GSL_2D_Spline ddD2LCDMddy_of_scale_spline;
#endif
#endif

  GSL_Spline D_spline;
  GSL_Spline dDdy_spline;   
  GSL_Spline ddDddy_spline; 
  GSL_Spline D2_spline;     
  GSL_Spline dD2dy_spline;  
  GSL_Spline ddD2ddy_spline;

  GSL_Spline DLCDM_spline;
  GSL_Spline dDLCDMdy_spline;   
  GSL_Spline ddDLCDMddy_spline; 
  GSL_Spline D2LCDM_spline;     
  GSL_Spline dD2LCDMdy_spline;  
  GSL_Spline ddD2LCDMddy_spline;

#ifdef MBETAMODEL
  GSL_Spline phi_of_a_spline;
#endif

} SplineContainer;

#if defined(MBETAMODEL)
//===============================================
// The integrand for phi(a) defined above
// Note dx = da/a so only a^3 in the denominator
//===============================================
double integrand_phiofa(double x, void *params){
  double a = exp(x);

  double beta  = beta_of_a(a);
  double mass2 = mass2_of_a(a);
  double source = 0.0;
  if(beta != 0.0){
    if(mass2 == 0.0){
      printf("Error in integrating phi(a) for MBETA models. m^2(a) = 0 and beta(a) != 0\n");
      MPI_Abort(MPI_COMM_WORLD,1);
      exit(1);
    }

    // This is dphi(a) / dlog(a)
    source = 9.0 * Omega * beta / (mass2 * a * a * a);
  }

  return source;
}

//===============================================
// Compute phi(a)
//===============================================
void compute_phi_of_a(){
  const int npts    = 1000;
  const double amin = 1.0/(1.0+REDSHIFT_START_INTEGRATION);
  const double amax = 1.0/(1.0+REDSHIFT_END_INTEGRATION);
  const double xmin = log(amin);
  const double xmax = log(amax);
  const double deltax = (xmax-xmin)/(double)(npts-1);
  const double phi_ini = 0.0;

  if(ThisTask == 0){
    printf("Compute screening theshold Phi_critical(a) from (m,beta)\n");
    printf("Assuming Phi_ini = %f at aini = %f\n", phi_ini, amin);
  }

  gsl_function F;
  F.function = &integrand_phiofa;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  double *x_arr   = my_malloc(npts * sizeof(double));
  double *phi_arr = my_malloc(npts * sizeof(double));
  double *err_arr = my_malloc(npts * sizeof(double));

  phi_arr[0] = phi_ini;
  x_arr[0] = xmin;
  err_arr[0] = 0.0;
  for(int i = 1; i < npts; i++){
    double xold = xmin + (i-1) * deltax;
    double xnow = xmin + i * deltax;

    double phi_local, error_local;
    gsl_integration_qag(&F, (double)xold, (double)xnow, 0, 1e-6, 1000, 6, w, &phi_local, &error_local); 

    x_arr[i] = xnow;
    phi_arr[i] = phi_arr[i-1] + phi_local;
    err_arr[i] = err_arr[i-1] + error_local;

  }

  // Rescale to phi_crit = phi(a)/2beta(a). For f(R) this becomes 1.5 f_R(a)
  for(int i = 0; i < npts; i++){
    double anow = exp(x_arr[i]);

    double beta = beta_of_a(anow);
    if(beta != 0.0){
      phi_arr[i] = phi_arr[i] / (2.0 * beta);
      err_arr[i] = err_arr[i] / (2.0 * beta);
    }

    if(ThisTask == 0 && i % 25 == 0){
      printf("a = %6.3f    Phi_critical = %10.5e   [Est. error in integration = %10.5e %%]\n", anow, phi_arr[i], (err_arr[i] / (phi_arr[i] + 1e-100)*100.0));

    }
  }

  // Make the spline
  Create_GSL_Spline(&SplineContainer.phi_of_a_spline, x_arr, phi_arr, npts);

  // Free up memory
  gsl_integration_workspace_free (w);
  my_free(x_arr);
  my_free(phi_arr);
  my_free(err_arr);
}

//================================================
// Spline lookup of the function computed above
//================================================
double phi_of_a_from_spline(double a){
  return Lookup_GSL_Spline(&SplineContainer.phi_of_a_spline, log(a));
}

#endif

double Qfactor(double a) { 
  return hubble(a)*a*a*a;
}

//==================================================================================
// First and second order growth equation. Params contains the k-value if it's 
// scale-dependent and a spline of the first order growth-factor
//==================================================================================
int ode_growth_D(double x, const double D[], double dDdy[], void *params){
  double a  = exp(x);
  double k  = *(double *) params;
  double mu = GeffoverG(a, k);
  double H  = hubble(a);
  double dH = dhubbleda(a);
  double betafac  = 1.5 * Omega * mu / ( a * a * a * H * H );
  double alphafac = 2.0 + a * dH / H;

#ifdef SCALEDEPENDENT
  // If we want to include the effects of gamma2 in this equation
  // Choose between [k, k/sq2, k/sq2, 0.0] or [k, k, 0.0, 0.0] or [k, k, k, 0.0]
  double modfac = (1.0 + 2.0 * pow2(a * a * H) * second_order_kernel(k, k/sqrt(2), k/sqrt(2), 0.0, a) / (1.5 * Omega * a * mu)); 
#else
  // If we have only time-dependent factors
  double modfac  = Factor_2LPT(a);
#endif
  
  //=======================================
  // Adding in effects of massive neutrinos
  //=======================================
  double munu_1LPT = GeffoverG_neutrino_1LPT(a, k);
  double munu_2LPT = GeffoverG_neutrino_2LPT(a, k);

  // First order growth factor
  dDdy[0] = D[1];
  dDdy[1] = - alphafac * D[1] + betafac * munu_1LPT * D[0];

  // Second order growth factor
  dDdy[2] = D[3];
  dDdy[3] = - alphafac * D[3] + betafac * munu_2LPT * (D[2] - modfac * D[0] * D[0] );

  return GSL_SUCCESS;
}

//==================================================================================
// First and second order LCDM growth equation. Params contains the k-value if it's 
// scale-dependent and a spline of the first order growth-factor
//==================================================================================
int ode_growth_DLCDM(double x, const double D[], double dDdy[], void *params){
  double a  = exp(x);
  double H  = hubble(a);
  double dH = dhubbleda(a);
  double k  = *(double *) params;
  double betafac = 1.5 * Omega / ( a * a * a * H * H );
  double alphafac = 2.0 + a * dH / H;

  //=======================================
  // Adding in effects of massive neutrinos
  //=======================================
  double munu_1LPT = GeffoverG_neutrino_1LPT(a, k);
  double munu_2LPT = GeffoverG_neutrino_2LPT(a, k);

  // First order growth factor
  dDdy[0] = D[1];
  dDdy[1] = - alphafac * D[1] + betafac * munu_1LPT * D[0];

  // Second order growth factor
  dDdy[2] = D[3];
  dDdy[3] = - alphafac * D[3] + betafac * munu_2LPT *  (D[2] - D[0] * D[0] );

  return GSL_SUCCESS;
}

#ifdef SCALEDEPENDENT

//==================================================================================
// Second order growth equation for the kernel k*L. Params contains the k,k1,k2-values 
// if it's scale-dependent and a spline of the first order growth-factor
//==================================================================================

int ode_second_order_growth_kernel_D2(double x, const double D2[], double dD2dy[], void *params){
  struct ode_second_order_growth_parameters *param_struct = (struct ode_second_order_growth_parameters*) params;
  double a     = exp(x);
  double k     = param_struct->k_value;
  double k1    = param_struct->k1_value;
  double k2    = param_struct->k2_value; // For comparing to approx use [k2 = k1 = k] for testing
  double cost  = param_struct->costheta_value;
  double H     = hubble(a);
  double dH    = dhubbleda(a);
  double mu_k  = GeffoverG(a, k);
  double mu_k1 = GeffoverG(a, k1);
  double mu_k2 = GeffoverG(a, k2);
  double betafac    = 1.5 * Omega / ( a * a * a * H * H );
  double alphafac   = 2.0 + a * dH / H;
  double factor2LPT = 1.0 - (mu_k1 + mu_k2 - mu_k)/mu_k * cost * cost + 2.0 * second_order_kernel(k, k1, k2, cost, a) / (betafac * mu_k);
  
  //=======================================
  // Adding in effects of massive neutrinos
  //=======================================
  double munu_1LPT = GeffoverG_neutrino_1LPT(a, k);
  double munu_2LPT = GeffoverG_neutrino_2LPT(a, k);

  // Second order equation
  dD2dy[0] = D2[1];
  dD2dy[1] = - alphafac * D2[1] + betafac * mu_k  * munu_2LPT * (D2[0] - factor2LPT * D2[2] * D2[4] );

  // First order equation for k1
  dD2dy[2] = D2[3];
  dD2dy[3] = - alphafac * D2[3] + betafac * mu_k1 * munu_1LPT * D2[2];

  // First order equation for k2
  dD2dy[4] = D2[5];
  dD2dy[5] = - alphafac * D2[5] + betafac * mu_k2 * munu_1LPT * D2[4];

  return GSL_SUCCESS;
}
#endif

//====================================================
// First order growth-factor and LCDM fitting-formula
//====================================================
double growth_D(double a){
  return Lookup_GSL_Spline(&SplineContainer.D_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_DLCDM(double a){
  return Lookup_GSL_Spline(&SplineContainer.DLCDM_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_D_LCDMFit(double a){
  return growthD(a);
}

//==================================================================
// Derivative of first order growth-factor and LCDM fitting-formula
//==================================================================
double growth_dDdy(double a){
  return Lookup_GSL_Spline(&SplineContainer.dDdy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D_spline,    log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_dDLCDMdy(double a){
  return Lookup_GSL_Spline(&SplineContainer.dDLCDMdy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,    log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_dDdy_LCDMFit(double a){
  return DprimeQ(a); 
}

//=========================================================================
// Second derivative of first order growth-factor and LCDM fitting-formula
//=========================================================================
double growth_ddDddy(double a){
  return Lookup_GSL_Spline(&SplineContainer.ddDddy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_ddDLCDMddy(double a){
  return Lookup_GSL_Spline(&SplineContainer.ddDLCDMddy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_ddDddy_LCDMFit(double a){
  return 1.5 * Omega * growthD(a) * a;
}

//====================================================
// Second order growth-factor and LCDM fitting-formula
//====================================================
double growth_D2(double a){
  return Lookup_GSL_Spline(&SplineContainer.D2_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D2_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_D2LCDM(double a){
  return Lookup_GSL_Spline(&SplineContainer.D2LCDM_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D2LCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_D2_LCDMFit(double a){
  return growthD2(a);
}

//==================================================================
// Derivative of second order growth-factor and LCDM fitting-formula
//==================================================================
double growth_dD2dy(double a){
  return Lookup_GSL_Spline(&SplineContainer.dD2dy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D2_spline,    log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_dD2LCDMdy(double a){
  return Lookup_GSL_Spline(&SplineContainer.dD2LCDMdy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D2LCDM_spline,    log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_dD2dy_LCDMFit(double a){
  return growthD2v(a); 
}

//===============================================================================
// Second order derivative of second order growth-factor and LCDM fitting-formula
//===============================================================================
double growth_ddD2ddy(double a){
  return Lookup_GSL_Spline(&SplineContainer.ddD2ddy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D2_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_ddD2LCDMddy(double a){
  return Lookup_GSL_Spline(&SplineContainer.ddD2LCDMddy_spline, log(a)) 
       / Lookup_GSL_Spline(&SplineContainer.D2LCDM_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR));
}
double growth_ddD2ddy_LCDMFit(double a){
  return 1.5 * Omega * pow2(growthD(a)) * (1.0 + 7./3.* pow( Omega / (Omega + (1 - Omega) * a * a * a), 1./143.)) * a;
}

//==================================================================
// This routine integrates the ODEs for the growth factors
// and splines it up
//==================================================================
void solve_for_growth_factors(){
  const int    npts   = 1000;
  const double zini   = REDSHIFT_START_INTEGRATION;
  const double zend   = REDSHIFT_END_INTEGRATION;
  const double aini   = 1.0/(1.0 + zini);
  const double xini   = log(aini);
  const double xend   = log(1.0/(1.0 + zend));
  const double deltax = (xend - xini) / (double) (npts-1);

  const int debug_show_comparison_to_fitting_func = 0;

  // Allocate arrays for growth-factors
  double *x_arr            = my_malloc(sizeof(double) * npts);

  // Modified gravity arrays
  double *D_arr            = my_malloc(sizeof(double) * npts);
  double *dDdy_arr         = my_malloc(sizeof(double) * npts);
  double *ddDddy_arr       = my_malloc(sizeof(double) * npts);
  double *D2_arr           = my_malloc(sizeof(double) * npts);
  double *dD2dy_arr        = my_malloc(sizeof(double) * npts);
  double *ddD2ddy_arr      = my_malloc(sizeof(double) * npts);

  // LCDM arrays
  double *DLCDM_arr        = my_malloc(sizeof(double) * npts);
  double *dDLCDMdy_arr     = my_malloc(sizeof(double) * npts);
  double *ddDLCDMddy_arr   = my_malloc(sizeof(double) * npts);
  double *D2LCDM_arr       = my_malloc(sizeof(double) * npts);
  double *dD2LCDMdy_arr    = my_malloc(sizeof(double) * npts);
  double *ddD2LCDMddy_arr  = my_malloc(sizeof(double) * npts);

  //=======================================================================
  // First order equation we solve below:
  // D''(x) + D'(x) ( 2  + Exp[x] H'[x]/H[x] ) - 1.5 * Omega * GeffG * exp(-3x) / H^2[x] D(x)
  // Defining the variable q = D'[x] then we get the system coupled first order system
  // dq[x]/dx = 1.5 * GeffG(x) * exp(-3x) / H^2[x] D[x] - q[x] ( 2  + Exp[x] H'[x]/H[x] )
  // dD[x]/dx = q[x]
  //=======================================================================

  //=======================================================================
  // Second order equation we need to solve (for LCDM)
  // D2''[x] + (2.0 + H'[x]/H[x]) D2'[x] = 3/2 1/a^3 1/H[x]^2 (D2[x] - D1[x]^2)
  // 
  // In the EdS we approx have D2 ~ -3/7 D1^2 and D1 ~ a in the matter era so
  // D2_ini ~ -3/7 and dD2_ini/dx ~ 2*D2_ini
  //
  // NB: Since we normalize at z=0 the factor -3/7 needs to be multiplied 
  // in to the initial displacement-field when used.
  //=======================================================================

  // Set up ODE system
  double k_value = 1e-3;
  gsl_odeiv2_system sys_D       = {ode_growth_D,     NULL, 4, &k_value};
  gsl_odeiv2_system sys_DLCDM   = {ode_growth_DLCDM, NULL, 4, &k_value};
  gsl_odeiv2_driver * ode_D     = gsl_odeiv2_driver_alloc_y_new (&sys_D,     gsl_odeiv2_step_rk2, MY_GSL_HSTART, MY_GSL_EPS, MY_GSL_REL);
  gsl_odeiv2_driver * ode_DLCDM = gsl_odeiv2_driver_alloc_y_new (&sys_DLCDM, gsl_odeiv2_step_rk2, MY_GSL_HSTART, MY_GSL_EPS, MY_GSL_REL);

  // Initial conditions for growing mode D ~ a so q ~ a
  double ode_D_x      = xini;
  double ode_DLCDM_x  = xini;
  double D_now[4]     = { 1.0, 1.0, -3.0/7.0, -6.0/7.0 };
  double DLCDM_now[4] = { 1.0, 1.0, -3.0/7.0, -6.0/7.0 };
  double muini        = GeffoverG(aini, k_value);

  // Store the initial values in array
  x_arr[0]           = xini;

  D_arr[0]           = D_now[0];
  dDdy_arr[0]        = D_now[1] * Qfactor(aini) / aini;
  ddDddy_arr[0]      = 1.5 * muini * Omega * aini * D_now[0];

  D2_arr[0]          = D_now[2];
  dD2dy_arr[0]       = D_now[3] * Qfactor(aini) / aini;
  ddD2ddy_arr[0]     = 1.5 * muini * Omega * aini * (D_now[2] - pow2( D_now[0] ) * Factor_2LPT(aini) );

  DLCDM_arr[0]       = DLCDM_now[0];
  dDLCDMdy_arr[0]    = DLCDM_now[1] * Qfactor(aini) / aini;
  ddDLCDMddy_arr[0]  = 1.5 * Omega * aini * DLCDM_now[0];

  D2LCDM_arr[0]      = DLCDM_now[2];
  dD2LCDMdy_arr[0]   = DLCDM_now[3] * Qfactor(aini) / aini;
  ddD2LCDMddy_arr[0] = 1.5 * Omega * aini * (DLCDM_now[2] - pow2( DLCDM_now[0] ) );

  // Integration in time
  for(int i = 1; i < npts; i++){
    double xnow = xini + i * deltax;
    double anow = exp(xnow);
    double mu   = GeffoverG(anow, k_value);

    // Integrate up MG growthfactor
    int status = gsl_odeiv2_driver_apply(ode_D, &ode_D_x, xnow, D_now);
    if(status != GSL_SUCCESS){
      printf("Error in integrating first order growth factor at x = %f  D = %f\n", xnow, D_now[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    // Integrate up LCDM growth factor
    int status_LCDM = gsl_odeiv2_driver_apply(ode_DLCDM, &ode_DLCDM_x, xnow, DLCDM_now);
    if(status_LCDM != GSL_SUCCESS){
      printf("Error in integrating first order growth factor for LCDM at x = %f  D = %f\n", xnow, DLCDM_now[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    // Store values
    x_arr[i]           = xnow;

    D_arr[i]           = D_now[0];
    dDdy_arr[i]        = D_now[1] * Qfactor(anow) / anow;
    ddDddy_arr[i]      = 1.5 * mu * Omega * anow * D_now[0];

    D2_arr[i]          = D_now[2];
    dD2dy_arr[i]       = D_now[3] * Qfactor(anow) / anow;
    ddD2ddy_arr[i]     = 1.5 * mu * Omega * anow * (D_now[2] - pow2( D_now[0] ) * Factor_2LPT(anow) );

    // Store values for LCDM
    DLCDM_arr[i]       = DLCDM_now[0];
    dDLCDMdy_arr[i]    = DLCDM_now[1] * Qfactor(anow) / anow;
    ddDLCDMddy_arr[i]  = 1.5 * Omega * anow * DLCDM_now[0];

    D2LCDM_arr[i]      = DLCDM_now[2];
    dD2LCDMdy_arr[i]   = DLCDM_now[3] * Qfactor(anow) / anow;
    ddD2LCDMddy_arr[i] = 1.5 * Omega * anow * (DLCDM_now[2] - pow2( DLCDM_now[0] ));
  }

  Create_GSL_Spline(&SplineContainer.D_spline,           x_arr, D_arr,           npts);
  Create_GSL_Spline(&SplineContainer.dDdy_spline,        x_arr, dDdy_arr,        npts);
  Create_GSL_Spline(&SplineContainer.ddDddy_spline,      x_arr, ddDddy_arr,      npts);
  Create_GSL_Spline(&SplineContainer.D2_spline,          x_arr, D2_arr,          npts);
  Create_GSL_Spline(&SplineContainer.dD2dy_spline,       x_arr, dD2dy_arr,       npts);
  Create_GSL_Spline(&SplineContainer.ddD2ddy_spline,     x_arr, ddD2ddy_arr,     npts);

  Create_GSL_Spline(&SplineContainer.DLCDM_spline,       x_arr, DLCDM_arr,       npts);
  Create_GSL_Spline(&SplineContainer.dDLCDMdy_spline,    x_arr, dDLCDMdy_arr,    npts);
  Create_GSL_Spline(&SplineContainer.ddDLCDMddy_spline,  x_arr, ddDLCDMddy_arr,  npts);
  Create_GSL_Spline(&SplineContainer.D2LCDM_spline,      x_arr, D2LCDM_arr,      npts);
  Create_GSL_Spline(&SplineContainer.dD2LCDMdy_spline,   x_arr, dD2LCDMdy_arr,   npts);
  Create_GSL_Spline(&SplineContainer.ddD2LCDMddy_spline, x_arr, ddD2LCDMddy_arr, npts);

  // Testing the new routines to the old routines (should match almost perfectly when Omega = 1.0)
  if(ThisTask == 0 && debug_show_comparison_to_fitting_func){
    int nout = 10;

    // NB: LCDM here only means no GeffG term so if hubble(a), dhubbleda(a) is modified then these will not match
    
    printf("\n==================================================================================================\n");
    printf("First order growth-factor for LCDM and derivatives compared to fitting functions used originally: \n");
    printf("==================================================================================================\n");
    for(int i = 0; i < nout; i++){
      double xnow = xini + (xend - xini)*i/(double)(nout-1);
      double anow = exp(xnow);
      printf(" a = %10.5f    D [%10.5f  %10.5f  Diff: %6.2f%%]    dD [%10.5f   %10.5f  Error: %6.2f%%]    ddD [%10.5f  %10.5f  Diff: %6.2f%%]\n", anow, 
          growth_DLCDM(anow),      growth_D_LCDMFit(anow),      fabs(growth_DLCDM(anow)/growth_D_LCDMFit(anow)-1.0)*100., 
          growth_dDLCDMdy(anow),   growth_dDdy_LCDMFit(anow),   fabs(growth_dDLCDMdy(anow)/growth_dDdy_LCDMFit(anow)-1.0)*100.,
          growth_ddDLCDMddy(anow), growth_ddDddy_LCDMFit(anow), fabs(growth_ddDLCDMddy(anow)/growth_ddDddy_LCDMFit(anow)-1.0)*100.);
    }

    printf("\n==================================================================================================\n");
    printf("Second order growth-factor for LCDM and derivatives compared to fitting functions used originally: \n");
    printf("==================================================================================================\n");
    for(int i = 0; i < nout; i++){
      double xnow = xini + (xend - xini)*i/(double)(nout-1);
      double anow = exp(xnow);
      printf(" a = %10.5f    D2[%10.5f  %10.5f  Diff: %6.2f%%]    dD2[%10.5f   %10.5f  Error: %6.2f%%]    ddD2[%10.5f  %10.5f  Diff: %6.2f%%]\n", anow, 
          growth_D2LCDM(anow),      growth_D2_LCDMFit(anow),      fabs(growth_D2LCDM(anow)/growth_D2_LCDMFit(anow)-1.0)*100., 
          growth_dD2LCDMdy(anow),   growth_dD2dy_LCDMFit(anow),   fabs(growth_dD2LCDMdy(anow)/growth_dD2dy_LCDMFit(anow)-1.0)*100.,
          growth_ddD2LCDMddy(anow), growth_ddD2ddy_LCDMFit(anow), fabs(growth_ddD2LCDMddy(anow)/growth_ddD2ddy_LCDMFit(anow)-1.0)*100.);
    }
  }

  // Print the modified gravity enhancement of the growth-factor
  if(ThisTask == 0 && modified_gravity_active && !use_lcdm_growth_factors){
    int nout = 10;
    printf("\n==================================================================================================\n");
    printf("First order growth-factor in MG wrt LCDM: \n");
    printf("==================================================================================================\n");
    for(int i = 0; i < nout; i++){
      double xnow = xini + (xend - xini)*i/(double)(nout-1);
      double anow = exp(xnow);
      printf(" a = %10.5f    D [%10.5f  %10.5f  Err: %6.2f%%]    dD [%10.5f   %10.5f  Error: %6.2f%%]    ddD [%10.5f  %10.5f  Err: %6.2f%%]\n", anow, 
          growth_D(anow),          growth_D_LCDMFit(anow),        fabs(growth_D(anow)/growth_D_LCDMFit(anow)-1.0)*100., 
          growth_dDdy(anow),       growth_dDdy_LCDMFit(anow),     fabs(growth_dDdy(anow)/growth_dDdy_LCDMFit(anow)-1.0)*100.,
          growth_ddDddy(anow),     growth_ddDddy_LCDMFit(anow),   fabs(growth_ddDddy(anow)/growth_ddDddy_LCDMFit(anow)-1.0)*100.);
    }

    printf("\n==================================================================================================\n");
    printf("Second order growth-factor in MG wrt LCDM: \n");
    printf("==================================================================================================\n");
    for(int i = 0; i < nout; i++){
      double xnow = xini + (xend - xini)*i/(double)(nout-1);
      double anow = exp(xnow);
      printf(" a = %10.5f    D2[%10.5f  %10.5f  Err: %6.2f%%]    dD2[%10.5f   %10.5f  Error: %6.2f%%]    ddD2[%10.5f  %10.5f  Err: %6.2f%%]\n", anow, 
          growth_D2(anow),         growth_D2_LCDMFit(anow),       fabs(growth_D2(anow)/growth_D2_LCDMFit(anow)-1.0)*100., 
          growth_dD2dy(anow),      growth_dD2dy_LCDMFit(anow),    fabs(growth_dD2dy(anow)/growth_dD2dy_LCDMFit(anow)-1.0)*100.,
          growth_ddD2ddy(anow),    growth_ddD2ddy_LCDMFit(anow),  fabs(growth_ddD2ddy(anow)/growth_ddD2ddy_LCDMFit(anow)-1.0)*100.);
    }
  }

  // Free up modified gravity arrays
  my_free(D_arr);
  my_free(dDdy_arr);
  my_free(ddDddy_arr);
  my_free(D2_arr);
  my_free(dD2dy_arr);
  my_free(ddD2ddy_arr);

  // Free up LCDM arrays
  my_free(DLCDM_arr);
  my_free(dDLCDMdy_arr);
  my_free(ddDLCDMddy_arr);
  my_free(D2LCDM_arr);
  my_free(dD2LCDMdy_arr);
  my_free(ddD2LCDMddy_arr);

  // Free up x-array
  my_free(x_arr);

#ifdef SCALEDEPENDENT

  // If we have scale-dependent growth-factors compute it here
  calculate_scale_dependent_growth_factor();

  //==================================================
  // To check the accuracy of our approximation
  // for the second order growth-factor
  //==================================================
  // 
  // check_error_approx();
  // MPI_Abort(MPI_COMM_WORLD,1); exit(1);
  //
  //==================================================

#endif

}

// Free up splines
void free_up_splines(){

  Free_GSL_Spline(&SplineContainer.D_spline);
  Free_GSL_Spline(&SplineContainer.dDdy_spline);
  Free_GSL_Spline(&SplineContainer.ddDddy_spline);
  Free_GSL_Spline(&SplineContainer.D2_spline);
  Free_GSL_Spline(&SplineContainer.dD2dy_spline);
  Free_GSL_Spline(&SplineContainer.ddD2ddy_spline);

  Free_GSL_Spline(&SplineContainer.DLCDM_spline);
  Free_GSL_Spline(&SplineContainer.dDLCDMdy_spline);
  Free_GSL_Spline(&SplineContainer.ddDLCDMddy_spline);
  Free_GSL_Spline(&SplineContainer.D2LCDM_spline);
  Free_GSL_Spline(&SplineContainer.dD2LCDMdy_spline);
  Free_GSL_Spline(&SplineContainer.ddD2LCDMddy_spline);

#ifdef SCALEDEPENDENT
  Free_GSL_2D_Spline(&SplineContainer.D_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.dDdy_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.ddDddy_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.dD2dy_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.ddD2ddy_of_scale_spline);
#ifdef MASSIVE_NEUTRINOS
  Free_GSL_2D_Spline(&SplineContainer.DLCDM_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.dDLCDMdy_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.ddDLCDMddy_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.D2LCDM_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.dD2LCDMdy_of_scale_spline);
  Free_GSL_2D_Spline(&SplineContainer.ddD2LCDMddy_of_scale_spline);
#endif
#endif

#ifdef MBETAMODEL
  Free_GSL_Spline(&SplineContainer.phi_of_a_spline);
#endif

}

//=================================================================
// The linear power-spectrum ratio P(k,a) / P_LCDM(k,a)
//=================================================================
double mg_pofk_ratio(double k, double a){
  if(! modified_gravity_active ) return 1.0;

#ifdef SCALEDEPENDENT
#ifdef MASSIVE_NEUTRINOS
  return pow2( Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,     log(a), log(k))
             / Lookup_GSL_2D_Spline(&SplineContainer.DLCDM_of_scale_spline, log(a), log(k)) );
#else
  return pow2( Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline, log(a), log(k))
             / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,         log(a)        ) );
#endif
#else
  return pow2( Lookup_GSL_Spline(&SplineContainer.D_spline,             log(a)        )        
             / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,         log(a)        ) );
#endif

}

//=================================================================
// Compute sigma8 enhancement at z=0 relative to LCDM
// This is useful when generating IC for MG starting from a LCDM
// P(k,z=0) such that we have exactly the same IC as for LCDM
// at the initial redshift (run with input_pofk_is_for_lcdm=1
// and input_sigma8_is_for_lcdm=1)
//=================================================================
double mg_sigma8_enhancement(double a){
  if(! modified_gravity_active ) return 1.0;

#ifndef SCALEDEPENDENT

  // If growth-factor is only time-dependent then no need to integrate
  return Lookup_GSL_Spline(&SplineContainer.D_spline,     log(a)) 
      /  Lookup_GSL_Spline(&SplineContainer.DLCDM_spline, log(a));

#else

  //=================================================================
  // sigma8 = Int P(k) * growth_factor^2 * k^3 dlog(k) * W^2(know*R8)
  //=================================================================

  // In units of h/Mpc
  const double kmin = 0.001;
  const double kmax = 100.0;
  const int npts    = 1000;

  double integrand      = 0.0;
  double integrand_LCDM = 0.0;
  double dlogk = log(kmax/kmin)/ (double)(npts-1);;
  for(int i = 0; i < npts; i++){
    double know  = exp( log(kmin) + log(kmax/kmin) * i / (double)(npts-1) );
    double kR8   = know * 8.0;
    double w     = 3.0/(kR8*kR8*kR8) * (sin(kR8) - kR8 * cos(kR8));
    double D     = Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline, log(a), log(know));
#ifdef MASSIVE_NEUTRINOS
    double DLCDM = Lookup_GSL_2D_Spline(&SplineContainer.DLCDM_of_scale_spline, log(a), log(know));
#else
    double DLCDM = Lookup_GSL_Spline(&SplineContainer.DLCDM_spline, log(a));
#endif
    double power = PowerSpec(know);

    integrand      += power * w * w * (D / DLCDM) * (D / DLCDM)   * know * know * know * dlogk / (2.0 * M_PI * M_PI);
    integrand_LCDM += power * w * w * know * know * know * dlogk / (2.0 * M_PI * M_PI);
  }
  integrand = sqrt(integrand);
  integrand_LCDM = sqrt(integrand_LCDM);
  return integrand / integrand_LCDM;

#endif
}

#ifdef SCALEDEPENDENT

//==============================================================================================
// Integrate the second order scale-dependent growth ODE for a given k-value (in units of h/Mpc)
//==============================================================================================
void integrate_scale_dependent_growth_ode(double *x_arr, double *d_arr, double *q_arr, double *d2_arr, double *q2_arr, double know, int npts){

  // Set up ODE system
  gsl_odeiv2_system sys_D       = {ode_growth_D, NULL, 4, &know};
  gsl_odeiv2_driver * ode_D     = gsl_odeiv2_driver_alloc_y_new (&sys_D, gsl_odeiv2_step_rk2, MY_GSL_HSTART, MY_GSL_EPS, MY_GSL_REL);

  // Set initial conditions (assuming they have been set in d[0] and q[0])
  double ode_D_x  = x_arr[0];
  double D_now[4] = { d_arr[0], q_arr[0], d2_arr[0], q2_arr[0] };

  // Integrate
  for(int i = 1; i < npts; i++){

    int status = gsl_odeiv2_driver_apply(ode_D, &ode_D_x, x_arr[i], D_now);
    if(status != GSL_SUCCESS){
      printf("Error in integrating scale-dependent growth factor at x = %f  D = %f\n", x_arr[i], D_now[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

    // Store values
    d_arr[i]  = D_now[0];
    q_arr[i]  = D_now[1];
    d2_arr[i] = D_now[2];
    q2_arr[i] = D_now[3];
  }

  // Free up memory
  gsl_odeiv2_driver_free(ode_D);
}

//==============================================================
// Computes the scale-dependent growth-factor for all values
// of k and stores it in splines
//==============================================================
void calculate_scale_dependent_growth_factor(){
  const int    npts = 1000;
  const double zini = REDSHIFT_START_INTEGRATION;
  const double zend = REDSHIFT_END_INTEGRATION;
  const double aini = 1.0/(1.0 + zini);
  const double xini = log(aini);
  const double xend = log(1.0/(1.0 + zend));

  // Limits for k-space growth-factor. k is here in units of h/Mpc
  const double kmin = 2.0 * M_PI / Box * 0.5;
  const double kmax = 2.0 * M_PI / Box * Nmesh * sqrt(3.0) * 2.0;
  const int    nk   = 1000;

  // Allocate memory
  double *logk_arr    = my_malloc(sizeof(double) * nk);
  double *x_arr       = my_malloc(sizeof(double) * npts);
  double *D_k_x       = my_malloc(sizeof(double) * npts * nk);
  double *dDdy_k_x    = my_malloc(sizeof(double) * npts * nk);
  double *ddDddy_k_x  = my_malloc(sizeof(double) * npts * nk);
  double *D2_k_x      = my_malloc(sizeof(double) * npts * nk);
  double *dD2dy_k_x   = my_malloc(sizeof(double) * npts * nk);
  double *ddD2ddy_k_x = my_malloc(sizeof(double) * npts * nk);

  // Define logk-array
  for(int i = 0; i < nk; i++)
    logk_arr[i] = log(kmin) + log(kmax/kmin) * i / (double) (nk-1);

  // Define x-array to store values in
  for(int i = 0; i < npts; i++)
    x_arr[i] = xini + (xend - xini) * i/(double) (npts-1);

#ifdef MASSIVE_NEUTRINOS

  //=======================================================
  // For massive neutrinos we have scale-dependent growth
  // for LCDM so...
  // Solve the growth ODEs for each value of k for LCDM
  //=======================================================

  int tmp_modified_gravity_active = modified_gravity_active;
  modified_gravity_active = 0;

  double *DLCDM_k_x       = malloc(sizeof(double) * npts * nk);
  double *dDLCDMdy_k_x    = malloc(sizeof(double) * npts * nk);
  double *ddDLCDMddy_k_x  = malloc(sizeof(double) * npts * nk);
  double *D2LCDM_k_x      = malloc(sizeof(double) * npts * nk);
  double *dD2LCDMdy_k_x   = malloc(sizeof(double) * npts * nk);
  double *ddD2LCDMddy_k_x = malloc(sizeof(double) * npts * nk);

  for(int k = 0; k < nk; k++){
    // Current value of k in h/Mpc
    double know = exp( log(kmin) + log(kmax/kmin) * k / (double) (nk-1) );

    double *D_arr       = &DLCDM_k_x[k * npts];
    double *dDdy_arr    = &dDLCDMdy_k_x[k * npts];
    double *ddDddy_arr  = &ddDLCDMddy_k_x[k * npts];
    double *D2_arr      = &D2LCDM_k_x[k * npts];
    double *dD2dy_arr   = &dD2LCDMdy_k_x[k * npts];
    double *ddD2ddy_arr = &ddD2LCDMddy_k_x[k * npts];

    // IC for growing mode
    D_arr[0]     = 1.0;
    dDdy_arr[0]  = 1.0;
    D2_arr[0]    = -3.0/7.0;
    dD2dy_arr[0] = -6.0/7.0;

    // Integrate PDE
    integrate_scale_dependent_growth_ode(x_arr, D_arr, dDdy_arr, D2_arr, dD2dy_arr, know, npts);

    // Set dDdy and ddDddy
    for(int i = 0; i < npts; i++){
      double anow    = exp(x_arr[i]);
      double mu      = GeffoverG(anow, know);
      dDdy_arr[i]    = dDdy_arr[i] * Qfactor(anow) / anow;
      ddDddy_arr[i]  = 1.5 * mu * Omega * anow * D_arr[i];
      dD2dy_arr[i]   = dD2dy_arr[i] * Qfactor(anow) / anow;
      ddD2ddy_arr[i] = 1.5 * mu * Omega * anow * (D2_arr[i] - pow2(D_arr[i]) );
    }
  }

  Create_GSL_2D_Spline(&SplineContainer.DLCDM_of_scale_spline,       x_arr, logk_arr, DLCDM_k_x,       npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.dDLCDMdy_of_scale_spline,    x_arr, logk_arr, dDLCDMdy_k_x,    npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.ddDLCDMddy_of_scale_spline,  x_arr, logk_arr, ddDLCDMddy_k_x,  npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.D2LCDM_of_scale_spline,      x_arr, logk_arr, D2LCDM_k_x,      npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.dD2LCDMdy_of_scale_spline,   x_arr, logk_arr, dD2LCDMdy_k_x,   npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.ddD2LCDMddy_of_scale_spline, x_arr, logk_arr, ddD2LCDMddy_k_x, npts, nk);

  free(DLCDM_k_x);
  free(dDLCDMdy_k_x);
  free(ddDLCDMddy_k_x);
  free(D2LCDM_k_x);
  free(dD2LCDMdy_k_x);
  free(ddD2LCDMddy_k_x);

  modified_gravity_active = tmp_modified_gravity_active;

#endif

  // Solve the growth ODEs for each value of k
  for(int k = 0; k < nk; k++){
    // Current value of k in h/Mpc
    double know = exp( log(kmin) + log(kmax/kmin) * k / (double) (nk-1) );

    double *D_arr       = &D_k_x[k * npts];
    double *dDdy_arr    = &dDdy_k_x[k * npts];
    double *ddDddy_arr  = &ddDddy_k_x[k * npts];
    double *D2_arr      = &D2_k_x[k * npts];
    double *dD2dy_arr   = &dD2dy_k_x[k * npts];
    double *ddD2ddy_arr = &ddD2ddy_k_x[k * npts];

    // IC for growing mode
    D_arr[0]     = 1.0;
    dDdy_arr[0]  = 1.0;
    D2_arr[0]    = -3.0/7.0;
    dD2dy_arr[0] = -6.0/7.0;

    // Integrate PDE
    integrate_scale_dependent_growth_ode(x_arr, D_arr, dDdy_arr, D2_arr, dD2dy_arr, know, npts);

    // Set dDdy and ddDddy
    for(int i = 0; i < npts; i++){
      double anow    = exp(x_arr[i]);
      double mu      = GeffoverG(anow, know);
      double munu_1LPT = GeffoverG_neutrino_1LPT(anow, know);
      double munu_2LPT = GeffoverG_neutrino_2LPT(anow, know);
      dDdy_arr[i]    = dDdy_arr[i] * Qfactor(anow) / anow;
      ddDddy_arr[i]  = 1.5 * mu * munu_1LPT * Omega * anow * D_arr[i];
      dD2dy_arr[i]   = dD2dy_arr[i] * Qfactor(anow) / anow;
      ddD2ddy_arr[i] = 1.5 * mu * munu_2LPT * Omega * anow * (D2_arr[i] - pow2(D_arr[i]) );
    }
  }

  Create_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,       x_arr, logk_arr, D_k_x,       npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.dDdy_of_scale_spline,    x_arr, logk_arr, dDdy_k_x,    npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.ddDddy_of_scale_spline,  x_arr, logk_arr, ddDddy_k_x,  npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline,      x_arr, logk_arr, D2_k_x,      npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.dD2dy_of_scale_spline,   x_arr, logk_arr, dD2dy_k_x,   npts, nk);
  Create_GSL_2D_Spline(&SplineContainer.ddD2ddy_of_scale_spline, x_arr, logk_arr, ddD2ddy_k_x, npts, nk);

  // Output some info about the splines
  if(ThisTask == 0 && modified_gravity_active && !use_lcdm_growth_factors){
    printf("\n===============================================\n");
    printf("Scale-dependent growth-factor relative to LCDM: \n");
    printf("===============================================\n");
    double anow = 1.0;
    for(int k = 0; k < nk; k+=100){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double dnow  = Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,  log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,          log(anow));
      double d2now = Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline, log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.D2LCDM_spline,         log(anow));
      double ratio = - 7.0 / 3.0 * Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline, log(anow), log(know))
                           / pow2( Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,  log(anow), log(know)) );
      printf("k = %8.3f h/Mpc  a = 1 at D / DLCDM           = %8.3f   D2 / D2LCDM           = %8.3f   D2 / ( - 3/7 D1^2 ) = %8.3f\n", know, dnow, d2now, ratio);
    }
    for(int k = 0; k < nk; k+=100){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double  dnow = Lookup_GSL_2D_Spline(&SplineContainer.dDdy_of_scale_spline,  log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.dDLCDMdy_spline,          log(anow));
      double d2now = Lookup_GSL_2D_Spline(&SplineContainer.dD2dy_of_scale_spline, log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.dD2LCDMdy_spline,         log(anow));
      printf("k = %8.3f h/Mpc  a = 1 at dDdy / dDLCDMdy     = %8.3f   dD2dy / dD2LCDMdy     = %8.3f\n", know, dnow, d2now);
    }
    for(int k = 0; k < nk; k+=100){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double  dnow = Lookup_GSL_2D_Spline(&SplineContainer.ddDddy_of_scale_spline,  log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.ddDLCDMddy_spline,          log(anow));
      double d2now = Lookup_GSL_2D_Spline(&SplineContainer.ddD2ddy_of_scale_spline, log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.ddD2LCDMddy_spline,         log(anow));
      printf("k = %8.3f h/Mpc  a = 1 at ddDddy / ddDLCDMddy = %8.3f   ddD2ddy / ddD2LCDMddy = %8.3f\n", know, dnow, d2now);
    }
    printf("\nSigma8(z=0) is enhanced by a factor: %f\n\n", mg_sigma8_enhancement(anow));
  }

  // Free up memory
  my_free(x_arr);
  my_free(logk_arr);
  my_free(D_k_x);  
  my_free(dDdy_k_x);
  my_free(ddDddy_k_x);
  my_free(D2_k_x);
  my_free(dD2dy_k_x);
  my_free(ddD2ddy_k_x);
}

//==================================================================================
// Scale-dependent growth-factor. Normalized to unity at a = 1 for all modes
// Assumes k in units of h/Mpc
//==================================================================================
double growth_D_scaledependent(double k, double a){
  return Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,       log(a),                                 log(k))
       / Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,       log(GROWTH_NORMALIZATION_SCALE_FACTOR), log(k));
}
double growth_dDdy_scaledependent(double k, double a){
  return Lookup_GSL_2D_Spline(&SplineContainer.dDdy_of_scale_spline,    log(a),                                 log(k))
       / Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,       log(GROWTH_NORMALIZATION_SCALE_FACTOR), log(k));
}
double growth_ddDddy_scaledependent(double k, double a){
  return Lookup_GSL_2D_Spline(&SplineContainer.ddDddy_of_scale_spline,  log(a),                                 log(k))
       / Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,       log(GROWTH_NORMALIZATION_SCALE_FACTOR), log(k));
}
double growth_D2_scaledependent(double k, double a){
  return Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline,      log(a),                                 log(k))
       / Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR), log(k));
}
double growth_dD2dy_scaledependent(double k, double a){
  return Lookup_GSL_2D_Spline(&SplineContainer.dD2dy_of_scale_spline,   log(a),                                 log(k))
       / Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR), log(k));
}
double growth_ddD2ddy_scaledependent(double k, double a){
  return Lookup_GSL_2D_Spline(&SplineContainer.ddD2ddy_of_scale_spline, log(a),                                 log(k))
       / Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR), log(k));
}

//==================================================================
// For computing single values of the growth kernel at [a]
// assuming D2 = 1 at a = aini 
// Returns the ratio wrt the approximate equation
//==================================================================

double compute_single_value_second_order_growth_kernel(double k, double k1, double k2, double a, int verbose){
  const int    npts = 1000;
  const double zini = REDSHIFT_START_INTEGRATION;
  const double aini = 1.0/(1.0 + zini);
  const double xini = log(aini);
  const double xend = log(a);
  const double deltax = (xend - xini)/(double)(npts-1);

  double cosphi = 1.0 - (k2*k2 - (k1-k)*(k1-k))/(2.0*k*k1);
  double costheta = -(k1 - k*cosphi)/k2;

  struct ode_second_order_growth_parameters ode_D2_kernel_param;
  ode_D2_kernel_param.k_value   = k;
  ode_D2_kernel_param.k1_value  = k1;
  ode_D2_kernel_param.k2_value  = k2;
  ode_D2_kernel_param.costheta_value = costheta;

  gsl_odeiv2_system sys_D2_kernel = {ode_second_order_growth_kernel_D2, NULL, 6, &ode_D2_kernel_param};
  gsl_odeiv2_driver * ode_D2_kernel = gsl_odeiv2_driver_alloc_y_new (&sys_D2_kernel, gsl_odeiv2_step_rk2, MY_GSL_HSTART, MY_GSL_EPS, MY_GSL_REL);

  // Initial conditions
  double ode_D2_x  = xini;
  double D2_now[6] = { -3.0/7.0 *(1.0 - costheta*costheta), -3.0/7.0 * 2.0 * (1.0 - costheta*costheta), 1.0, 1.0, 1.0, 1.0 };

  // Now we can integrate over k
  for(int i = 1; i < npts; i++){
    double xnow = xini + i * deltax;

    // Integrate up MG growthfactor
    int status = gsl_odeiv2_driver_apply(ode_D2_kernel, &ode_D2_x, xnow, D2_now);
    if(status != GSL_SUCCESS){
      printf("Error in integrating second order growth kernel at x = %f  D2 = %f\n", xnow, D2_now[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
  }

  if(ThisTask == 0 && verbose)
    printf("k: %8.3f  k1: %8.3f  k2: %8.3f  costheta: %8.3f  D1: %8.3f  D2: %8.3f\n", k, k1, k2, costheta, D2_now[2], D2_now[0]);

  double D2_approx = Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline, log(a), log(k)) * (1.0 - costheta*costheta); 
  double D2_exact  = D2_now[0];

  return D2_approx / D2_exact;
}

//==================================================================
// Compute the error in the approximative equation for D2 for
// different triangle configurations
//==================================================================
void check_error_approx(){
  double a = 1.0;
  double k,k1,k2,res;
  int verbose = 0;

  // kmin/kmax are in units of h/Mpc
  const double kmin = 2.0 * M_PI / Box * 0.5;
  const double kmax = 2.0 * M_PI / Box * Nmesh * sqrt(3.0) * 2.0;
  int npts = 30;

  printf("Equilateral: k = k1 = k2\n");
  for(int i = 0; i < npts; i++){

    k = k1 = k2 = exp(log(kmin) + log(kmax/kmin)*i/(double)(npts-1));
    res = compute_single_value_second_order_growth_kernel(k,k1,k2,a,verbose);
    printf("%f    %f\n", k , res);
  }
  
  printf("\nSqueezed k = k1  and k2 = 0\n");
  for(int i = 0; i < npts; i++){

    k = k1 = exp(log(kmin) + log(kmax/kmin)*i/(double)(npts-1));
    k2 = 0.001;
    res = compute_single_value_second_order_growth_kernel(k,k1,k2,a,verbose);
    printf("%f    %f\n", k , res);
  }
  
  printf("\nOrthogonal k1 = k2 = k/sqrt(2)\n");
  for(int i = 0; i < npts; i++){

    k = exp(log(kmin) + log(kmax/kmin)*i/(double)(npts-1));
    k1 = k2 = k/sqrt(2.0);
    res = compute_single_value_second_order_growth_kernel(k,k1,k2,a,verbose);
    printf("%f    %f\n", k , res);
  }
}
#endif

//=============================================================
// Below follow the original PICOLA functions
//=============================================================

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
// 
// Deprecated methods from the original PICOLA code
// These have been replaced. Only used to compare to the new version
//
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

