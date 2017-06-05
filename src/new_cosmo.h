#ifndef NEWCOSMOINC
#define NEWCOSMOINC

//==========================================================================//
//                                                                          //
//  MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017           //
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
  Create_GSL_Spline(&SplineContainer.phi_of_a_spline, x_arr, phi_arr, npts);

  // Free up memory
  gsl_integration_workspace_free (w);
  free(x_arr);
  free(phi_arr);
  free(err_arr);
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
  // check_error_approx();
  // exit(1);
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
  return pow2( Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline, log(a), log(k)) 
             / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,         log(a)        ) );
#else
  return pow2( Lookup_GSL_Spline(&SplineContainer.D_spline,             log(a)        )        
             / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,         log(a)        ) );
#endif

}

//=================================================================
// Compute sigma8 enhancement at z=0 relative to LCDM
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
    double DLCDM = Lookup_GSL_Spline(&SplineContainer.DLCDM_spline, log(a));
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
  double *logk_arr    = malloc(sizeof(double) * nk);
  double *x_arr       = malloc(sizeof(double) * npts);
  double *D_k_x       = malloc(sizeof(double) * npts * nk);
  double *dDdy_k_x    = malloc(sizeof(double) * npts * nk);
  double *ddDddy_k_x  = malloc(sizeof(double) * npts * nk);
  double *D2_k_x      = malloc(sizeof(double) * npts * nk);
  double *dD2dy_k_x   = malloc(sizeof(double) * npts * nk);
  double *ddD2ddy_k_x = malloc(sizeof(double) * npts * nk);

  // Define logk-array
  for(int i = 0; i < nk; i++)
    logk_arr[i] = log(kmin) + log(kmax/kmin) * i / (double) (nk-1);

  // Define x-array to store values in
  for(int i = 0; i < npts; i++)
    x_arr[i] = xini + (xend - xini) * i/(double) (npts-1);

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
  if(ThisTask == 0 && modified_gravity_active){
    printf("\n===============================================\n");
    printf("Scale-dependent growth-factor relative to LCDM: \n");
    printf("===============================================\n");
    double anow = 1.0;
    for(int k = 0; k < nk; k+=50){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double dnow  = Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,  log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.DLCDM_spline,          log(anow));
      double d2now = Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline, log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.D2LCDM_spline,         log(anow));
      double ratio = - 7.0 / 3.0 * Lookup_GSL_2D_Spline(&SplineContainer.D2_of_scale_spline, log(anow), log(know))
                           / pow2( Lookup_GSL_2D_Spline(&SplineContainer.D_of_scale_spline,  log(anow), log(know)) );
      printf("k = %8.3f h/Mpc  a = 1 at D / DLCDM           = %8.3f   D2 / D2LCDM           = %8.3f   D2 / ( - 3/7 D1^2 ) = %8.3f\n", know, dnow, d2now, ratio);
    }
    for(int k = 0; k < nk; k+=50){
      double know = exp( log(kmin) + log(kmax/kmin) * k/(double)(nk-1) );
      double  dnow = Lookup_GSL_2D_Spline(&SplineContainer.dDdy_of_scale_spline,  log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.dDLCDMdy_spline,          log(anow));
      double d2now = Lookup_GSL_2D_Spline(&SplineContainer.dD2dy_of_scale_spline, log(anow), log(know)) 
                   / Lookup_GSL_Spline(&SplineContainer.dD2LCDMdy_spline,         log(anow));
      printf("k = %8.3f h/Mpc  a = 1 at dDdy / dDLCDMdy     = %8.3f   dD2dy / dD2LCDMdy     = %8.3f\n", know, dnow, d2now);
    }
    for(int k = 0; k < nk; k+=50){
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
  free(x_arr);
  free(logk_arr);
  free(D_k_x);  
  free(dDdy_k_x);
  free(ddDddy_k_x);
  free(D2_k_x);
  free(dD2dy_k_x);
  free(ddD2ddy_k_x);
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

#endif
