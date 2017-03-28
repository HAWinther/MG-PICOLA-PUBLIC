//==========================================================================//
//                                                                          //
// MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017            //
//                                                                          //
// This file contains updated cosmology functions needed to solve the       //
// general growth ODEs plus look-up functions                               //
//                                                                          //
// The original functions used in PICOLA are found as                       //
// [ xxxx_LCDMFit(double a) ]                                               //
//                                                                          //
// Below d/dy = Q(a) d/da                                                   //
//                                                                          //
//==========================================================================//

// Some useful functions
#define MIN(x,y) ((x) < (y) ? (y) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define pow2(x)  ((x)*(x))
#define pow3(x)  ((x)*(x)*(x))
#define pow4(x)  ((x)*(x)*(x)*(x))
#define pow5(x)  ((x)*(x)*(x)*(x)*(x))

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

#include "Spline.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include "user_defined_functions.h"

#if defined(MBETAMODEL)
//===============================================
// The integrand for phi(a) defined above
//===============================================
double integrand_phiofa(double x, void *params){
  double a = exp(x);
  
  // This is dphi(a) / dlog(a)
  double source = 9.0 * Omega * beta_of_a(a) / (mass2_of_a(a) * a * a * a);
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

  double *x_arr   = malloc(npts * sizeof(double));
  double *phi_arr = malloc(npts * sizeof(double));
  double *err_arr = malloc(npts * sizeof(double));

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

    phi_arr[i] = phi_arr[i] / (2.0 * beta_of_a(anow));
    err_arr[i] = err_arr[i] / (2.0 * beta_of_a(anow));

    if(ThisTask == 0 && i % 100 == 0){
      printf("a = %6.3f    Phi_critical = %10.5e   [Est. error in integration = %10.5e %%]\n", anow, phi_arr[i], (err_arr[i] / phi_arr[i])*100.0);

    }
  }

  // Make the spline
  phi_of_a_spline = malloc(sizeof(Spline));
  create_spline(phi_of_a_spline,  x_arr,  phi_arr, npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);

  // Free up memory
  gsl_integration_workspace_free (w);
  free(x_arr);
  free(phi_arr);
  free(err_arr);
}
#endif

double Qfactor(double a) { 
  return hubble(a)*a*a*a;
}

//========================================================
// Splines needed to store growth-factors and derivatives
//========================================================

#ifdef SCALEDEPENDENT
// Scale-dependent growth-factors
SplineArray FirstOrderGrowthFactor_D;
SplineArray FirstOrderGrowthFactor_dDdy;
SplineArray FirstOrderGrowthFactor_ddDddy;

SplineArray SecondOrderGrowthFactor_D;
SplineArray SecondOrderGrowthFactor_dDdy;
SplineArray SecondOrderGrowthFactor_ddDddy;
#endif

struct TimeDependentSplineContainer{
  int Splines_Are_Created;

  Spline *D_spline;
  Spline *dDdy_spline;
  Spline *ddDddy_spline;
  Spline *D2_spline;
  Spline *dD2dy_spline;
  Spline *ddD2ddy_spline;

  Spline *DLCDM_spline;
  Spline *dDLCDMdy_spline;
  Spline *ddDLCDMddy_spline;
  Spline *D2LCDM_spline;
  Spline *dD2LCDMdy_spline;
  Spline *ddD2LCDMddy_spline;

} TimeDependentSplines;

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

  // First order growth factor
  dDdy[0] = D[1];
  dDdy[1] = - alphafac * D[1] + betafac * D[0];

  // Second order growth factor
  dDdy[2] = D[3];
  dDdy[3] = - alphafac * D[3] + betafac * (D[2] - modfac * D[0] * D[0] );
  
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
  double betafac = 1.5 * Omega / ( a * a * a * H * H );
  double alphafac = 2.0 + a * dH / H;

  // First order growth factor
  dDdy[0] = D[1];
  dDdy[1] = - alphafac * D[1] + betafac * D[0];

  // Second order growth factor
  dDdy[2] = D[3];
  dDdy[3] = - alphafac * D[3] + betafac * (D[2] - D[0] * D[0] );

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
  double factor2LPT = 1.0 - cost * cost + 2.0 * second_order_kernel(k, k1, k2, cost, a) / (betafac * mu_k);

  // Second order equation
  dD2dy[0] = D2[1];
  dD2dy[1] = - alphafac * D2[1] + betafac * mu_k  * (D2[0] - factor2LPT * D2[2] * D2[4] );
  
  // First order equation for k1
  dD2dy[2] = D2[3];
  dD2dy[3] = - alphafac * D2[3] + betafac * mu_k1 * D2[2];
  
  // First order equation for k2
  dD2dy[4] = D2[5];
  dD2dy[5] = - alphafac * D2[5] + betafac * mu_k2 * D2[4];
  
  return GSL_SUCCESS;
}
#endif

//====================================================
// First order growth-factor and LCDM fitting-formula
//====================================================
double growth_D(double a){
  return spline_lookup(TimeDependentSplines.D_spline,           log(a) ) / spline_lookup(TimeDependentSplines.D_spline,     log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_DLCDM(double a){
  return spline_lookup(TimeDependentSplines.DLCDM_spline,       log(a) ) / spline_lookup(TimeDependentSplines.DLCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_D_LCDMFit(double a){
  return growthD(a);
}

//==================================================================
// Derivative of first order growth-factor and LCDM fitting-formula
//==================================================================
double growth_dDdy(double a){
  return spline_lookup(TimeDependentSplines.dDdy_spline,        log(a) ) / spline_lookup(TimeDependentSplines.D_spline,     log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_dDLCDMdy(double a){
  return spline_lookup(TimeDependentSplines.dDLCDMdy_spline,    log(a) ) / spline_lookup(TimeDependentSplines.DLCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_dDdy_LCDMFit(double a){
  return DprimeQ(a); 
}

//=========================================================================
// Second derivative of first order growth-factor and LCDM fitting-formula
//=========================================================================
double growth_ddDddy(double a){
  return spline_lookup(TimeDependentSplines.ddDddy_spline,      log(a) ) / spline_lookup(TimeDependentSplines.D_spline,      log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_ddDLCDMddy(double a){
  return spline_lookup(TimeDependentSplines.ddDLCDMddy_spline,  log(a) ) / spline_lookup(TimeDependentSplines.DLCDM_spline,  log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_ddDddy_LCDMFit(double a){
  return 1.5 * Omega * growthD(a) * a;
}

//====================================================
// Second order growth-factor and LCDM fitting-formula
//====================================================
double growth_D2(double a){
  return spline_lookup(TimeDependentSplines.D2_spline,          log(a) ) / spline_lookup(TimeDependentSplines.D2_spline,     log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_D2LCDM(double a){
  return spline_lookup(TimeDependentSplines.D2LCDM_spline,      log(a) ) / spline_lookup(TimeDependentSplines.D2LCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_D2_LCDMFit(double a){
  return growthD2(a);
}

//==================================================================
// Derivative of second order growth-factor and LCDM fitting-formula
//==================================================================
double growth_dD2dy(double a){
  return spline_lookup(TimeDependentSplines.dD2dy_spline,       log(a) ) / spline_lookup(TimeDependentSplines.D2_spline,     log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_dD2LCDMdy(double a){
  return spline_lookup(TimeDependentSplines.dD2LCDMdy_spline,   log(a) ) / spline_lookup(TimeDependentSplines.D2LCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_dD2dy_LCDMFit(double a){
  return growthD2v(a); 
}

//===============================================================================
// Second order derivative of second order growth-factor and LCDM fitting-formula
//===============================================================================
double growth_ddD2ddy(double a){
  return spline_lookup(TimeDependentSplines.ddD2ddy_spline,     log(a) ) / spline_lookup(TimeDependentSplines.D2_spline,     log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
}
double growth_ddD2LCDMddy(double a){
  return spline_lookup(TimeDependentSplines.ddD2LCDMddy_spline, log(a) ) / spline_lookup(TimeDependentSplines.D2LCDM_spline, log(GROWTH_NORMALIZATION_SCALE_FACTOR) );
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

  // Allocate arrays for growth-factors
  double *x_arr            = malloc(sizeof(double) * npts);
  
  // Modified gravity arrays
  double *D_arr            = malloc(sizeof(double) * npts);
  double *dDdy_arr         = malloc(sizeof(double) * npts);
  double *ddDddy_arr       = malloc(sizeof(double) * npts);
  double *D2_arr           = malloc(sizeof(double) * npts);
  double *dD2dy_arr        = malloc(sizeof(double) * npts);
  double *ddD2ddy_arr      = malloc(sizeof(double) * npts);
  
  // LCDM arrays
  double *DLCDM_arr        = malloc(sizeof(double) * npts);
  double *dDLCDMdy_arr     = malloc(sizeof(double) * npts);
  double *ddDLCDMddy_arr   = malloc(sizeof(double) * npts);
  double *D2LCDM_arr       = malloc(sizeof(double) * npts);
  double *dD2LCDMdy_arr    = malloc(sizeof(double) * npts);
  double *ddD2LCDMddy_arr  = malloc(sizeof(double) * npts);

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
  double k_value = 0.0;
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
 
  // Make first order growth factor splines for modifed gravity
  TimeDependentSplines.D_spline            = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.D_spline,           x_arr, D_arr,           npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.dDdy_spline         = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.dDdy_spline,        x_arr, dDdy_arr,        npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.ddDddy_spline       = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.ddDddy_spline,      x_arr, ddDddy_arr,      npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  
  // Make first order growth factor splines for LCDM
  TimeDependentSplines.DLCDM_spline        = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.DLCDM_spline,       x_arr, DLCDM_arr,       npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.dDLCDMdy_spline     = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.dDLCDMdy_spline,    x_arr, dDLCDMdy_arr,    npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.ddDLCDMddy_spline   = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.ddDLCDMddy_spline,  x_arr, ddDLCDMddy_arr,  npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);

  // Make second order growth factor splines for modifed gravity
  TimeDependentSplines.D2_spline           = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.D2_spline,          x_arr, D2_arr,          npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.dD2dy_spline        = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.dD2dy_spline,       x_arr, dD2dy_arr,       npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.ddD2ddy_spline      = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.ddD2ddy_spline,     x_arr, ddD2ddy_arr,     npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  
  // Make second order growth factor splines for LCDM
  TimeDependentSplines.D2LCDM_spline       = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.D2LCDM_spline,      x_arr, D2LCDM_arr,      npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.dD2LCDMdy_spline    = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.dD2LCDMdy_spline,   x_arr, dD2LCDMdy_arr,   npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  TimeDependentSplines.ddD2LCDMddy_spline  = malloc(sizeof(Spline));
  create_spline(TimeDependentSplines.ddD2LCDMddy_spline, x_arr, ddD2LCDMddy_arr, npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
 
  // Set flag to signal that we have allocated memory for the splines
  TimeDependentSplines.Splines_Are_Created = 1;
 
  // Testing the new routines to the old routines (should match almost perfectly when Omega = 1.0)
  if(ThisTask == 0){
    int nout = 10;
    printf("\n==================================================================================================\n");
    printf("First order growth-factor for LCDM and derivatives compared to fitting functions used originally: \n");
    printf("==================================================================================================\n");
    for(int i = 0; i < nout; i++){
      double xnow = xini + (xend - xini)*i/(double)(nout-1);
      double anow = exp(xnow);
      printf(" a = %10.5f    D [%10.5f  %10.5f  Err: %6.2f%%]    dD [%10.5f   %10.5f  Error: %6.2f%%]    ddD [%10.5f  %10.5f  Err: %6.2f%%]\n", anow, 
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
      printf(" a = %10.5f    D2[%10.5f  %10.5f  Err: %6.2f%%]    dD2[%10.5f   %10.5f  Error: %6.2f%%]    ddD2[%10.5f  %10.5f  Err: %6.2f%%]\n", anow, 
          growth_D2LCDM(anow),      growth_D2_LCDMFit(anow),      fabs(growth_D2LCDM(anow)/growth_D2_LCDMFit(anow)-1.0)*100., 
          growth_dD2LCDMdy(anow),   growth_dD2dy_LCDMFit(anow),   fabs(growth_dD2LCDMdy(anow)/growth_dD2dy_LCDMFit(anow)-1.0)*100.,
          growth_ddD2LCDMddy(anow), growth_ddD2ddy_LCDMFit(anow), fabs(growth_ddD2LCDMddy(anow)/growth_ddD2ddy_LCDMFit(anow)-1.0)*100.);
    }
  }

  // Print the modified gravity enhancement of the growth-factor
  if(ThisTask == 0 && modified_gravity_active){
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
  free(D_arr);
  free(dDdy_arr);
  free(ddDddy_arr);
  free(D2_arr);
  free(dD2dy_arr);
  free(ddD2ddy_arr);
  
  // Free up LCDM arrays
  free(DLCDM_arr);
  free(dDLCDMdy_arr);
  free(ddDLCDMddy_arr);
  free(D2LCDM_arr);
  free(dD2LCDMdy_arr);
  free(ddD2LCDMddy_arr);
 
  // Free up x-array
  free(x_arr);

#ifdef SCALEDEPENDENT
  
  // If we have scale-dependent growth-factors compute it here
  calculate_scale_dependent_growth_factor();

  //==================================================
  // To check the accuracy of our approximation
  // for the second order growth-factor
  //==================================================
  // 
  // solve_for_second_order_kernel();
  // exit(1);
  // 
  // Or faster...
  // check_error_approx();
  // exit(1);
  //
  //==================================================

#endif

}

// Free up splines
void free_up_splines(){
  if(TimeDependentSplines.Splines_Are_Created == 1){
    
    // Free modified gravity splines
    free_spline(TimeDependentSplines.D_spline);
    free_spline(TimeDependentSplines.dDdy_spline);
    free_spline(TimeDependentSplines.ddDddy_spline);
    free_spline(TimeDependentSplines.D2_spline);
    free_spline(TimeDependentSplines.dD2dy_spline);
    free_spline(TimeDependentSplines.ddD2ddy_spline);
    free(TimeDependentSplines.D_spline);
    free(TimeDependentSplines.dDdy_spline);
    free(TimeDependentSplines.ddDddy_spline);
    free(TimeDependentSplines.D2_spline);
    free(TimeDependentSplines.dD2dy_spline);
    free(TimeDependentSplines.ddD2ddy_spline);

    // Free LCDM splines
    free_spline(TimeDependentSplines.DLCDM_spline);
    free_spline(TimeDependentSplines.dDLCDMdy_spline);
    free_spline(TimeDependentSplines.ddDLCDMddy_spline);
    free_spline(TimeDependentSplines.D2LCDM_spline);
    free_spline(TimeDependentSplines.dD2LCDMdy_spline);
    free_spline(TimeDependentSplines.ddD2LCDMddy_spline);
    free(TimeDependentSplines.DLCDM_spline);
    free(TimeDependentSplines.dDLCDMdy_spline);
    free(TimeDependentSplines.ddDLCDMddy_spline);
    free(TimeDependentSplines.D2LCDM_spline);
    free(TimeDependentSplines.dD2LCDMdy_spline);
    free(TimeDependentSplines.ddD2LCDMddy_spline);
  }

#if defined(MBETAMODEL)
  if(phi_of_a_spline != NULL){
    free_spline(phi_of_a_spline);
    free(phi_of_a_spline);
    phi_of_a_spline = NULL;
  }
#endif

#ifdef SCALEDEPENDENT
  // Free up scale-dependent first order growth-factor
  free_SplineArray(&FirstOrderGrowthFactor_D);
  free_SplineArray(&FirstOrderGrowthFactor_dDdy);
  free_SplineArray(&FirstOrderGrowthFactor_ddDddy);

  // Free up scale-dependent second order growth-factor
  free_SplineArray(&SecondOrderGrowthFactor_D);
  free_SplineArray(&SecondOrderGrowthFactor_dDdy);
  free_SplineArray(&SecondOrderGrowthFactor_ddDddy);
#endif
}

//=================================================================
// The linear power-spectrum ratio P(k,a) / P_LCDM(k,a)
//=================================================================
double mg_pofk_ratio(double k, double a){
  if(! modified_gravity_active ) return 1.0;
#ifdef SCALEDEPENDENT
  return pow2( interpolate_from_splinearray(&FirstOrderGrowthFactor_D, k, log(a)) / spline_lookup(TimeDependentSplines.DLCDM_spline, log(a)) );
#else
  return pow2( spline_lookup(TimeDependentSplines.D_spline, log(a)) / spline_lookup(TimeDependentSplines.DLCDM_spline, log(a)) );
#endif
}

//=================================================================
// Compute sigma8 enhancement at z=0 relative to LCDM
//=================================================================
double mg_sigma8_enhancement(double a){

#ifndef SCALEDEPENDENT

  // If growth-factor is only time-dependent then no need to integrate
  return spline_lookup(TimeDependentSplines.D_spline, log(a)) / spline_lookup(TimeDependentSplines.DLCDM_spline, log(a));

#else

  //=================================================================
  // sigma8 = Int P(k) * growth_factor^2 * k^3 dlog(k) * W^2(know*R8)
  //=================================================================

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
    double D     = interpolate_from_splinearray(&FirstOrderGrowthFactor_D,  know, log(a));
    double DLCDM = spline_lookup(TimeDependentSplines.DLCDM_spline, log(a));
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
// of k and stores it in SplineArray FirstOrderGrowthFactor_X
// and SecondOrderGrowthFactor_X
//==============================================================
void calculate_scale_dependent_growth_factor(){
  const int    npts = 1000;
  const double zini = REDSHIFT_START_INTEGRATION;
  const double zend = REDSHIFT_END_INTEGRATION;
  const double aini = 1.0/(1.0 + zini);
  const double xini = log(aini);
  const double xend = log(1.0/(1.0 + zend));

  // Limits for k-space growth-factor. Boxsize is assume to be in units of h/Mpc
  const double kmin = 2.0 * M_PI / Box * 0.5;
  const double kmax = 2.0 * M_PI / Box * Nmesh * sqrt(3.0) * 2.0;
  const int    nk   = 1000;

  // Make first order growth-factors
  FirstOrderGrowthFactor_D.kmin = kmin;
  FirstOrderGrowthFactor_D.kmax = kmax;
  FirstOrderGrowthFactor_D.nk   = nk;
  FirstOrderGrowthFactor_D.splinearray = malloc(sizeof(Spline *) * nk);
  FirstOrderGrowthFactor_D.is_created = 1;
  
  FirstOrderGrowthFactor_dDdy.kmin = kmin;
  FirstOrderGrowthFactor_dDdy.kmax = kmax;
  FirstOrderGrowthFactor_dDdy.nk   = nk;
  FirstOrderGrowthFactor_dDdy.splinearray = malloc(sizeof(Spline *) * nk);
  FirstOrderGrowthFactor_dDdy.is_created = 1;
  
  FirstOrderGrowthFactor_ddDddy.kmin = kmin;
  FirstOrderGrowthFactor_ddDddy.kmax = kmax;
  FirstOrderGrowthFactor_ddDddy.nk   = nk;
  FirstOrderGrowthFactor_ddDddy.splinearray = malloc(sizeof(Spline *) * nk);
  FirstOrderGrowthFactor_ddDddy.is_created = 1;

  // Make second order growth-factors
  SecondOrderGrowthFactor_D.kmin = kmin;
  SecondOrderGrowthFactor_D.kmax = kmax;
  SecondOrderGrowthFactor_D.nk   = nk;
  SecondOrderGrowthFactor_D.splinearray = malloc(sizeof(Spline *) * nk);
  SecondOrderGrowthFactor_D.is_created = 1;

  SecondOrderGrowthFactor_dDdy.kmin = kmin;
  SecondOrderGrowthFactor_dDdy.kmax = kmax;
  SecondOrderGrowthFactor_dDdy.nk   = nk;
  SecondOrderGrowthFactor_dDdy.splinearray = malloc(sizeof(Spline *) * nk);
  SecondOrderGrowthFactor_dDdy.is_created = 1;
  
  SecondOrderGrowthFactor_ddDddy.kmin = kmin;
  SecondOrderGrowthFactor_ddDddy.kmax = kmax;
  SecondOrderGrowthFactor_ddDddy.nk   = nk;
  SecondOrderGrowthFactor_ddDddy.splinearray = malloc(sizeof(Spline *) * nk);
  SecondOrderGrowthFactor_ddDddy.is_created = 1;
  
  // Allocate memory
  double *x_arr      = malloc(sizeof(double) * npts);
  
  double *D_arr      = malloc(sizeof(double) * npts);
  double *dDdy_arr   = malloc(sizeof(double) * npts);
  double *ddDddy_arr = malloc(sizeof(double) * npts);
  double *q_arr = dDdy_arr;
  
  double *D2_arr      = malloc(sizeof(double) * npts);
  double *dD2dy_arr   = malloc(sizeof(double) * npts);
  double *ddD2ddy_arr = malloc(sizeof(double) * npts);
  double *q2_arr = dD2dy_arr;

  // Define x-array to store values in
  for(int i = 0; i < npts; i++)
    x_arr[i] = xini + (xend - xini) * i/(double) (npts-1);

  // Solve the growth ODEs for each value of k
  for(int k = 0; k < nk; k++){
    // Current value of k in h/Mpc
    double know = exp( log(kmin) + log(kmax/kmin) * k / (double) (nk-1) );
    
    // IC for growing mode
    D_arr[0]     = 1.0;
    dDdy_arr[0]  = 1.0;
    D2_arr[0]    = -3.0/7.0;
    dD2dy_arr[0] = -6.0/7.0;
    
    // Integrate PDE
    integrate_scale_dependent_growth_ode(x_arr, D_arr, q_arr, D2_arr, q2_arr, know, npts);
    
    // Set dDdy and ddDddy
    for(int i = 0; i < npts; i++){
      double anow    = exp(x_arr[i]);
      double mu      = GeffoverG(anow, know);
      dDdy_arr[i]    = q_arr[i] * Qfactor(anow) / anow;
      ddDddy_arr[i]  = 1.5 * mu * Omega * anow * D_arr[i];
      dD2dy_arr[i]   = q2_arr[i] * Qfactor(anow) / anow;
      ddD2ddy_arr[i] = 1.5 * mu * Omega * anow * (D2_arr[i] - pow2(D_arr[i]) );
    }

    // Create splines
    FirstOrderGrowthFactor_D.splinearray[k]      = malloc(sizeof(Spline));
    create_spline(FirstOrderGrowthFactor_D.splinearray[k],       x_arr, D_arr,       npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
    FirstOrderGrowthFactor_dDdy.splinearray[k]   = malloc(sizeof(Spline));
    create_spline(FirstOrderGrowthFactor_dDdy.splinearray[k],    x_arr, dDdy_arr,    npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
    FirstOrderGrowthFactor_ddDddy.splinearray[k] = malloc(sizeof(Spline));
    create_spline(FirstOrderGrowthFactor_ddDddy.splinearray[k],  x_arr, ddDddy_arr,  npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
    
    SecondOrderGrowthFactor_D.splinearray[k]      = malloc(sizeof(Spline));
    create_spline(SecondOrderGrowthFactor_D.splinearray[k],      x_arr, D2_arr,      npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
    SecondOrderGrowthFactor_dDdy.splinearray[k]   = malloc(sizeof(Spline));
    create_spline(SecondOrderGrowthFactor_dDdy.splinearray[k],   x_arr, dD2dy_arr,   npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
    SecondOrderGrowthFactor_ddDddy.splinearray[k] = malloc(sizeof(Spline));
    create_spline(SecondOrderGrowthFactor_ddDddy.splinearray[k], x_arr, ddD2ddy_arr, npts, BC_NATURAL_SPLINE, BC_NATURAL_SPLINE, LINEAR_SPACED_SPLINE);
  }

  // Output some info about the splines
  if(ThisTask == 0){
    printf("\n===============================================\n");
    printf("Scale-dependent growth-factor relative to LCDM: \n");
    printf("===============================================\n");
    double anow = 1.0;
    for(int k = 0; k < nk; k+=50){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double dnow  = interpolate_from_splinearray(&FirstOrderGrowthFactor_D,  know, log(anow)) / spline_lookup(TimeDependentSplines.DLCDM_spline,  log(anow) );
      double d2now = interpolate_from_splinearray(&SecondOrderGrowthFactor_D, know, log(anow)) / spline_lookup(TimeDependentSplines.D2LCDM_spline, log(anow) );
      double ratio = - 7.0 / 3.0 * interpolate_from_splinearray(&SecondOrderGrowthFactor_D, know, log(anow)) 
                           / pow2( interpolate_from_splinearray(&FirstOrderGrowthFactor_D,  know, log(anow))); 
      printf("k = %8.3f h/Mpc    D(k, a=1):    : MG / LCDM = %8.3f  D2      : MG / LCDM = %8.3f   D2 / ( - 3/7 D1^2 ) = %8.3f\n", know, dnow, d2now, ratio);
    }
    for(int k = 0; k < nk; k+=50){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double dnow  = interpolate_from_splinearray(&FirstOrderGrowthFactor_dDdy,  know, log(anow)) / spline_lookup(TimeDependentSplines.dDLCDMdy_spline,  log(anow) );
      double d2now = interpolate_from_splinearray(&SecondOrderGrowthFactor_dDdy, know, log(anow)) / spline_lookup(TimeDependentSplines.dD2LCDMdy_spline, log(anow) );
      printf("k = %8.3f h/Mpc    dDdy(k, a=1)  : MG / LCDM = %8.3f  dD2dy   : MG / LCDM = %8.3f\n", know, dnow, d2now);
    }
    for(int k = 0; k < nk; k+=50){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double dnow  = interpolate_from_splinearray(&FirstOrderGrowthFactor_ddDddy,  know, log(anow)) / spline_lookup(TimeDependentSplines.ddDLCDMddy_spline,  log(anow) );
      double d2now = interpolate_from_splinearray(&SecondOrderGrowthFactor_ddDddy, know, log(anow)) / spline_lookup(TimeDependentSplines.ddD2LCDMddy_spline, log(anow) );
      printf("k = %8.3f h/Mpc    ddDddy(k, a=1): MG / LCDM = %8.3f  ddD2ddy : MG / LCDM = %8.3f\n", know, dnow, d2now);
    }
    printf("\nSigma8(z=0) is enhanced by a factor: %f\n\n", mg_sigma8_enhancement(anow));
  }

  // Free up memory
  free(x_arr);

  free(D_arr);
  free(dDdy_arr);
  free(ddDddy_arr);
  
  free(D2_arr);
  free(dD2dy_arr);
  free(ddD2ddy_arr);
}

//==================================================================================
// Scale-dependent growth-factor. Normalized to unity at a = 1 for all modes
// Assumes k in units of h/Mpc
//==================================================================================
double growth_D_scaledependent(double k, double a){
  return (   interpolate_from_splinearray(&FirstOrderGrowthFactor_D,       k, log(a) ) 
           / interpolate_from_splinearray(&FirstOrderGrowthFactor_D,       k, log(GROWTH_NORMALIZATION_SCALE_FACTOR) ) );
}
double growth_dDdy_scaledependent(double k, double a){
  return (   interpolate_from_splinearray(&FirstOrderGrowthFactor_dDdy,    k, log(a) ) 
           / interpolate_from_splinearray(&FirstOrderGrowthFactor_D,       k, log(GROWTH_NORMALIZATION_SCALE_FACTOR) ) );
}
double growth_ddDddy_scaledependent(double k, double a){
  return (   interpolate_from_splinearray(&FirstOrderGrowthFactor_ddDddy,  k, log(a) ) 
           / interpolate_from_splinearray(&FirstOrderGrowthFactor_D,       k, log(GROWTH_NORMALIZATION_SCALE_FACTOR) ) );
}
double growth_D2_scaledependent(double k, double a){
  return (   interpolate_from_splinearray(&SecondOrderGrowthFactor_D,      k, log(a) ) 
           / interpolate_from_splinearray(&SecondOrderGrowthFactor_D,      k, log(GROWTH_NORMALIZATION_SCALE_FACTOR) ) );
}
double growth_dD2dy_scaledependent(double k, double a){
  return (   interpolate_from_splinearray(&SecondOrderGrowthFactor_dDdy,   k, log(a) ) 
           / interpolate_from_splinearray(&SecondOrderGrowthFactor_D,      k, log(GROWTH_NORMALIZATION_SCALE_FACTOR) ) );
}
double growth_ddD2ddy_scaledependent(double k, double a){
  return (   interpolate_from_splinearray(&SecondOrderGrowthFactor_ddDddy, k, log(a) ) 
           / interpolate_from_splinearray(&SecondOrderGrowthFactor_D,      k, log(GROWTH_NORMALIZATION_SCALE_FACTOR) ) );
}

//=================================================================================
// L(k1,k2,k) = L{ |k1|, |k2|, cos(theta) }
// We use nk points in the k-directions + ncos points on the circle
// This means we need to solve nk * nk * nphi ODEs
//=================================================================================
void solve_for_second_order_kernel(){
  const int nk      = 101;
  const int nphi    = 101;
  const int npts    = 100;
  const double kmin = 2.0 * M_PI / Box * 0.5;
  const double kmax = 2.0 * M_PI / Box * Nmesh * sqrt(3.0) * 2.0;
  const double cmin = -1.0;
  const double cmax = 1.0;
  const double zini = REDSHIFT_START_INTEGRATION;
  const double zend = REDSHIFT_END_INTEGRATION;
  const double aini = 1.0/(1.0 + zini);
  const double xini = log(aini);
  const double xend = log(1.0/(1.0 + zend));
  const double deltax = (xend - xini) / (double) (npts-1);

  // For testing
  int ipresent_time_low  = (int)(-xini/deltax);
  int ipresent_time_high = ipresent_time_low + 1;
  double x_present_time_low  = xini + deltax * ipresent_time_low; 
  double x_present_time_high  = xini + deltax * ipresent_time_high; 
  printf("z = 0 is at a = %f\n", exp((x_present_time_high + x_present_time_low)/2.0));

  // Allocate grid [not free'd as written now...]
  InterpolationGrid *g = malloc(sizeof(InterpolationGrid));
  g->ntot     = nk * nk * nphi;
  g->n[0]     = nk;
  g->n[1]     = nk;
  g->n[2]     = nphi;
  g->logamin  = xini;
  g->logamax  = xend;
  g->xlow[0]  = g->xlow[1] = log(kmin);
  g->xlow[2]  = cmin;
  g->xhigh[0] = g->xhigh[1] = log(kmax);
  g->xhigh[2] = cmax;
  g->ntime    = npts;
  g->D2       = malloc(sizeof(double *) * g->ntime);
  g->dD2dy    = malloc(sizeof(double *) * g->ntime);
  g->ddD2ddy  = malloc(sizeof(double *) * g->ntime);
  for(int i = 0; i < npts; i++){
    g->D2[i]      = malloc(sizeof(double) * g->ntot);
    g->dD2dy[i]   = malloc(sizeof(double) * g->ntot);
    g->ddD2ddy[i] = malloc(sizeof(double) * g->ntot);
  }

  // Temporary memory
  double *tmp_D2      = malloc(sizeof(double)*npts);
  double *tmp_dD2dy   = malloc(sizeof(double)*npts);
  double *tmp_ddD2ddy = malloc(sizeof(double)*npts);

  // Loop over all combinations of |k1|, |k2|, |k| and phi(k,k1) 
  for(int ik = 0; ik < nk; ik++){
    printf("Computing for i = %i / %i\n", ik, nk);
    double k = exp(log(kmin) + log(kmax/kmin) * ik / (double)(nk-1));
    for(int ik1 = 0; ik1 < nk; ik1++){
      double k1 = exp(log(kmin) + log(kmax/kmin) * ik1 / (double)(nk-1));

      for(int ic = 0; ic < nphi; ic++){
        double cosphi = cmin + (cmax - cmin) * ic /(double)(nphi-1);

        // The current value of k2
        double k2 = sqrt((k1-k)*(k1-k) + 2.0*k*k1*(1.0 - cosphi));

        // The current value of costheta
        double costheta;
        if(k2 > 0.0){
          costheta = (k1 - k*cosphi)/k2;
        } else {
          costheta = 0.0;
        }

        // Roundoff can give costheta slightly larger than 1 or less than -1 so just in case fix it
        if(costheta < -1.0) costheta = -1.0;
        if(costheta >  1.0) costheta =  1.0;

        // Initialize parameters struct
        struct ode_second_order_growth_parameters ode_D2_kernel_param;
        ode_D2_kernel_param.k_value   = k;
        ode_D2_kernel_param.k1_value  = k1;
        ode_D2_kernel_param.k2_value  = k2;
        ode_D2_kernel_param.costheta_value = costheta;

        // Set up ODE system
        gsl_odeiv2_system sys_D2_kernel = {ode_second_order_growth_kernel_D2, NULL, 6, &ode_D2_kernel_param};
        gsl_odeiv2_driver * ode_D2_kernel = gsl_odeiv2_driver_alloc_y_new (&sys_D2_kernel, gsl_odeiv2_step_rk2, MY_GSL_HSTART_LOWACC, MY_GSL_EPS_LOWACC, MY_GSL_REL_LOWACC);

        // Initial conditions
        double ode_D2_x      = xini;
        double D2_now[6]     = { -3.0/7.0 *(1 - costheta*costheta), -3.0/7.0 * 2.0 * (1 - costheta*costheta), 1.0, 1.0, 1.0, 1.0 };

        // Stored the initial values in array
        tmp_D2[0]      = D2_now[0];
        tmp_dD2dy[0]   = D2_now[1] * Qfactor(aini) / aini;
        tmp_ddD2ddy[0] = 0.0;

        // Now we can integrate over k
        for(int i = 1; i < npts; i++){
          double xnow = xini + i * deltax;
          double anow = exp(xnow);

          // Integrate up MG growthfactor
          int status = gsl_odeiv2_driver_apply(ode_D2_kernel, &ode_D2_x, xnow, D2_now);
          if(status != GSL_SUCCESS){
            printf("Error in integrating second order growth kernel at x = %f  D2 = %f\n", xnow, D2_now[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1);
          }

          // Store values
          tmp_D2[i]      = D2_now[0];
          tmp_dD2dy[i]   = D2_now[1] * Qfactor(anow) / anow;
          tmp_ddD2ddy[i] = 0.0;
        }
 
        // Store data in interpolation-grid
        int index = ik + g->n[0] * (ik1 + g->n[1] * ic);
        for(int i = 0; i < npts; i++){
          g->D2[i][index]      = tmp_D2[i];
          g->dD2dy[i][index]   = tmp_dD2dy[i];
          g->ddD2ddy[i][index] = tmp_ddD2ddy[i];
        }

        // For testing
        if(ik==ik1 && (ic == (nphi-1)/2 || ic == 5 || ic == 50) ){
          double xx = xini + deltax * (ipresent_time_low);
          printf("cos = %f  Deff = %f    (%f)\n", costheta, tmp_D2[ipresent_time_low], 
              interpolate_from_splinearray(&SecondOrderGrowthFactor_D, k, xx)*(1.0 - costheta*costheta));
        }

        // Free up memory
        gsl_odeiv2_driver_free(ode_D2_kernel);
      }
    }
  }

  // Free up memory
  free(tmp_D2);
  free(tmp_dD2dy);
  free(tmp_ddD2ddy);
}

//===========================================================================================================
// Only working for 1CPU. The index we send in is the grid-index, not physical wavenumber.
// This is normalized such that if computes the convolution as computed by squaring the real density field
// and FFTing it.
// This can be used for computing and comparing the second order displacment-field compared to approximation
//===========================================================================================================

void integrate_up_kernel(int ik_grid_x, int ik_grid_y, int ik_grid_z, complex_kind *deltak_grid, float_kind *result){
  complex_kind integral;
  double scale = 2.0 * M_PI / Box;

  // Integer k-vector
  int ik_vec_x = (ik_grid_x > Nmesh/2 ? ik_grid_x - Nmesh : ik_grid_x);
  int ik_vec_y = (ik_grid_y > Nmesh/2 ? ik_grid_y - Nmesh : ik_grid_y);
  int ik_vec_z = (ik_grid_z > Nmesh/2 ? ik_grid_z - Nmesh : ik_grid_z);

  // k-vector in units of h/Mpc
  double k_vec_x = ik_vec_x * scale;
  double k_vec_y = ik_vec_y * scale;
  double k_vec_z = ik_vec_z * scale;
  double k_mag = sqrt(k_vec_x*k_vec_x + k_vec_y*k_vec_y + k_vec_z*k_vec_z);

  (void)(k_mag);

  double volume = 1.0;
  integral[0] = integral[1] = 0.0;
  for (int i = 0; i < Nmesh; i++){
    for (int j = 0; j < Nmesh; j++){
      for (int k = 0; k < Nmesh; k++){

        // k1-index in grid
        int ik1_grid_x = i;
        int ik1_grid_y = j;
        int ik1_grid_z = k;

        // k2-index in grid
        int ik2_grid_x = mymod(ik_grid_x - i, Nmesh);
        int ik2_grid_y = mymod(ik_grid_y - j, Nmesh);
        int ik2_grid_z = mymod(ik_grid_z - k, Nmesh);

        unsigned int ind_ik1 = (ik1_grid_x*Nmesh + ik1_grid_y)*Nmesh + ik1_grid_z;
        unsigned int ind_ik2 = (ik2_grid_x*Nmesh + ik2_grid_y)*Nmesh + ik2_grid_z;

        double DeltaDelta_Re = deltak_grid[ind_ik1][0]*deltak_grid[ind_ik2][0] - deltak_grid[ind_ik1][1]*deltak_grid[ind_ik2][1];
        double DeltaDelta_Im = deltak_grid[ind_ik1][0]*deltak_grid[ind_ik2][1] + deltak_grid[ind_ik1][1]*deltak_grid[ind_ik2][0];

        integral[0] += DeltaDelta_Re;
        integral[1] += DeltaDelta_Im;
      }
    }
  }
  integral[0] *= pow3(1.0 / (double) Nmesh);
  integral[1] *= pow3(1.0 / (double) Nmesh);

  printf("Volume: %f\n", volume);

  result[0] = integral[0];
  result[1] = integral[1];
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
  double costheta = (k1 - k*cosphi)/k2;

  struct ode_second_order_growth_parameters ode_D2_kernel_param;
  ode_D2_kernel_param.k_value   = k;
  ode_D2_kernel_param.k1_value  = k1;
  ode_D2_kernel_param.k2_value  = k2;
  ode_D2_kernel_param.costheta_value = costheta;

  gsl_odeiv2_system sys_D2_kernel = {ode_second_order_growth_kernel_D2, NULL, 6, &ode_D2_kernel_param};
  gsl_odeiv2_driver * ode_D2_kernel = gsl_odeiv2_driver_alloc_y_new (&sys_D2_kernel, gsl_odeiv2_step_rk2, MY_GSL_HSTART, MY_GSL_EPS, MY_GSL_REL);

  // Initial conditions
  double ode_D2_x      = xini;
  double D2_now[6]     = { -3.0/7.0 *(1 - costheta*costheta), -3.0/7.0 * 2.0 * (1 - costheta*costheta), 1.0, 1.0, 1.0, 1.0 };

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

  double D2_approx = interpolate_from_splinearray(&SecondOrderGrowthFactor_D, k, log(a) ) * (1.0 - costheta*costheta); 
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
  
  const double kmin = 2.0 * M_PI / Box * 0.5;
  const double kmax = 2.0 * M_PI / Box * Nmesh * sqrt(3.0) * 2.0;
  int npts = 30;

  // Equilateral
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
