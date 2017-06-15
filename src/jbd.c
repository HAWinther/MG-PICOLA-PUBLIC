#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include "wrappers.h"

#define JBD_EPSILON_CONVERGENCE 1e-6
#define JBD_REDSHIFT_START      1e6
#define JBD_NPOINTS             10000

GSL_Spline JBD_GeffSpline;
GSL_Spline JBD_HSpline;
GSL_Spline JBD_dHdaSpline;

//=============================================================================
// 
// This module takes in physical cosmological parameters and the JBD
// constant and solves the background equation and computes the hubble constant
// today.
//
// Input:    OmegaMh2, OmegaRh2, OmegaVh2, w, npts
// Returns:  HubbleParameter. 
//
// After we are done solving then H(a), dHda(a) and Geff(a) is availiable
// externally from the functions listed below
//
// void   JBD_Solve_Background(double w, double omegamh2, double omegavh2, double omegarh2, double *h);
// double JBD_Hubble_of_a(double a);
// double JBD_dHubbleda_of_a(double a);
// double JBD_GeffG_of_a(double a);
//
//=============================================================================

//=============================================================================
// List of internal functions
//=============================================================================

double JBD_HubbleFunction(double x, double y, double dy);
double JBD_dHubbleFunctiondx(double x, double y, double dy);
double JBD_HubbleFunction_of_z(double x, double y, double z);
double JBD_HLCDM(double x);
double JBD_OmegaM(double x, double y, double dy);
double JBD_OmegaR(double x, double y, double dy);
double JBD_OmegaPhi(double x, double y, double dy);
double JBD_GeffOverG(double y);
int    JBD_ode(double x, const double y[], double dydx[], void *params);
void   JBD_solve_ode(double yi, double dyi, double *x_arr, double *y_arr, double *dy_arr);

//=============================================================================
// Cosmological and model parameters
//=============================================================================
struct Parameters {
  double Omegam0h2; // kappa^2 rhom0/[3 H^2 phi0]
  double Omegar0h2; // kappa^2 rhor0/[3 H^2 phi0]
  double Omegav0h2; // kappa^2 Lambda/[3 H^2 phi0]
  double w;         // The JBD parameters

  double zini;      // Starting redshift for integration
  double epsilon;   // Convergence criterion
  double dy0;       // Value of dy/dx(x=0)
  double yi;        // Initial condition yi(zini)

  int npts;         // Number of points to store
  double *x_arr;    // x  = log(a) array
  double *y_arr;    // y  = log(phi/phi0) array
  double *dy_arr;   // dy = dlog(phi)/dloga array

  int ThisTask;     // MPI parameter
} gg;

//=============================================================================
// Lookup-functions for use after we have done the calculation
//=============================================================================
double JBD_Hubble_of_a(double a){
  return Lookup_GSL_Spline(&JBD_HSpline, log(a));
}
double JBD_dHubbleda_of_a(double a){
  return Lookup_GSL_Spline(&JBD_dHdaSpline, log(a));
}
double JBD_GeffG_of_a(double a){
  return Lookup_GSL_Spline(&JBD_GeffSpline, log(a));
}

//=============================================================================
// The Hubble function in terms of x = log(a), y = log(phi/phi0) and dy = dy/dx
// We have H(a=1) = h = little h
//=============================================================================
double JBD_HubbleFunction(double x, double y, double dy){
  return sqrt( exp(-y) * (gg.Omegar0h2 * exp(-4*x) + gg.Omegam0h2 * exp(-3*x) + gg.Omegav0h2) / (1.0 + dy - gg.w/6 * dy * dy) );
}
double JBD_dHubbleFunctiondx(double x, double y, double dy){
  double H   = JBD_HubbleFunction(x, y, dy);
  double om  = gg.Omegam0h2 * exp(-y - 3.0*x) / (H*H);
  double ov  = gg.Omegav0h2 * exp(-y) / (H*H);
  double ora = gg.Omegar0h2 * exp(-y - 4.0*x) / (H*H);
  return H/2.0 * (-3.0 - ora - gg.w/2.0*dy*dy + dy - 3.0*(om + 4.0*ov)/(3.0 + 2.0*gg.w) + 3.0*ov);
}

//=============================================================================
// The Hubble function in terms of x = log(a), y = log(phi/phi0) and z = H a^3 e^y dy/dx
// We have H(a=1) = h = little h
//=============================================================================
double JBD_HubbleFunction_of_z(double x, double y, double z){
  return sqrt( exp(-y) * (gg.Omegar0h2 * exp(-4*x) + gg.Omegam0h2 * exp(-3*x) + gg.Omegav0h2)  + z*z/12.0 * (2.0*gg.w + 3.0) * exp(-2*y-6*x) ) - z/2.0 * exp(-y-3*x);
}

//=============================================================================
// Hubble function for LCDM with the same value of Omega_m0 and Omega_r0
// We have H(a=1) = h = little h
//=============================================================================
double JBD_HLCDM(double x){
  return sqrt(gg.Omegar0h2 * exp(-4*x) + gg.Omegam0h2 * exp(-3*x) + gg.Omegav0h2);
}

//=============================================================================
// Density parameters
//=============================================================================
double JBD_OmegaM(double x, double y, double dy){
  return gg.Omegam0h2 * exp(-3*x) / pow( JBD_HubbleFunction(x, y, dy) , 2);
}
double JBD_OmegaR(double x, double y, double dy){
  return gg.Omegar0h2 * exp(-4*x) / pow( JBD_HubbleFunction(x, y, dy) , 2);
}
double JBD_OmegaPhi(double x, double y, double dy){
  return 1.0 - JBD_OmegaM(x,y,dy) - JBD_OmegaR(x,y,dy);
}

//=============================================================================
// GeffOverG = phi0/phi = e^{-y}
//=============================================================================
double JBD_GeffOverG(double y){
  return exp(-y);
}

//=============================================================================
// The ODE for the JBD scalar field
//=============================================================================
int JBD_ode(double x, const double y[], double dydx[], void *params){
  double hubb = JBD_HubbleFunction_of_z(x, y[0], y[1]);

  dydx[0] = y[1] * exp(- y[0] - 3*x) / hubb;
  dydx[1] = 3.0*exp(3.0*x) / (2*gg.w+3) / hubb * (gg.Omegam0h2 * exp(-3*x) + 4.0 * gg.Omegav0h2);

  return GSL_SUCCESS;
}

//=============================================================================
// Solve the ODE for the initial condition (y, dy/dx) = (yi, dyi) 
// Stores the values found in the provided arrays
//=============================================================================
void JBD_solve_ode(double yi, double dyi, double *x_arr, double *y_arr, double *dy_arr){
  const double xini = log(1.0 / (1.0 + gg.zini));
  const double xend = log(1.0);
  const double deltax = (xend-xini)/(double)(gg.npts-1);

  // Set up ODE system
  gsl_odeiv2_system JBD_sys = {JBD_ode, NULL, 2, NULL};
  gsl_odeiv2_driver * JBD_ode_driver = gsl_odeiv2_driver_alloc_y_new (&JBD_sys,  gsl_odeiv2_step_rk2, 1e-8, 1e-8, 0.0);

  // Set IC
  double y[2] = {yi, JBD_HubbleFunction(xini, yi, dyi) * exp(3.0*xini + yi) * dyi };

  double ode_x = xini;

  // Set IC in array
  x_arr[0]  = xini;
  y_arr[0]  = y[0];
  dy_arr[0] = y[1];

  // Solve the ODE
  for(int i = 1; i < gg.npts; i++){
    double xnow = xini + i * deltax;

    int status = gsl_odeiv2_driver_apply(JBD_ode_driver, &ode_x, xnow, y);
    if(status != GSL_SUCCESS){
      printf("Error in integrating at x = %f  y = %f\n", xnow, y[0]);
      exit(1);
    }

    x_arr[i]  = xnow;
    y_arr[i]  = y[0];
    dy_arr[i] = y[1] / JBD_HubbleFunction_of_z(xnow, y[0], y[1]) * exp(-3.0*xnow - y[0]);
  }
}

//=============================================================================
// Find correct initial conditions using bisection
//=============================================================================
void JBD_find_correct_IC_using_bisection(){
  const int npts = gg.npts;
  double *x_arr  = gg.x_arr;
  double *y_arr  = gg.y_arr;
  double *dy_arr = gg.dy_arr;

  // Find the correct initial condition
  double philow  = 0.0;
  double phihigh = 1.0;
  double phinow, yi, dyi;

  int istep = 0;
  while(1){
    ++istep;

    // Current value for phii
    phinow = (philow+phihigh)/2.0;

    // Solve ODE
    yi  = gg.yi = log(phinow);
    dyi = 0.0;
    JBD_solve_ode(yi, dyi, x_arr, y_arr, dy_arr);

    // Check for convergence
    double phiphi0 = exp(y_arr[npts-1]);
    if( fabs(phiphi0 - 1.0) < gg.epsilon ) {
#ifdef JBDDEBUG
      printf("Convergence of solution found after %i iterations\n", istep);
#endif
      break;
    }

    // Bisection step
    if(phiphi0 < 1.0){
      philow = phinow;
    } else {
      phihigh = phinow;
    }
  }
}

void JBD_Solve_Background(double w, double omegamh2, double omegavh2, double omegarh2, double *h){
  
  //=========================================
  // Set the parameters needed by the solver
  //=========================================
  gg.w         = w;
  gg.Omegam0h2 = omegamh2;
  gg.Omegar0h2 = omegarh2;
  gg.Omegav0h2 = omegavh2;
  gg.npts      = JBD_NPOINTS;
  gg.zini      = JBD_REDSHIFT_START;
  gg.epsilon   = JBD_EPSILON_CONVERGENCE;
 
  gg.x_arr     = my_malloc(sizeof(double)*gg.npts);
  gg.y_arr     = my_malloc(sizeof(double)*gg.npts);
  gg.dy_arr    = my_malloc(sizeof(double)*gg.npts);
  
  // Find the correct IC (after this is done we have the correct solution in gg-arrays)
  JBD_find_correct_IC_using_bisection();

  // The Hubble parameter we find
  *h = JBD_HubbleFunction(0.0, gg.y_arr[gg.npts-1], gg.dy_arr[gg.npts-1]);
  
  // Make arrays for splining
  double *loga = my_malloc(sizeof(double)*gg.npts);
  double *H    = my_malloc(sizeof(double)*gg.npts);
  double *dHda = my_malloc(sizeof(double)*gg.npts);
  double *Geff = my_malloc(sizeof(double)*gg.npts);
  for(int i = 0; i < gg.npts; i++){
    loga[i] = gg.x_arr[i];
    Geff[i] = JBD_GeffOverG(gg.y_arr[i]);
    H[i]    = JBD_HubbleFunction(gg.x_arr[i], gg.y_arr[i], gg.dy_arr[i]) / (*h);
    dHda[i] = JBD_dHubbleFunctiondx(gg.x_arr[i], gg.y_arr[i], gg.dy_arr[i])/(*h) * exp(-gg.x_arr[i]);
  }
  
  // Spline up results
  Create_GSL_Spline(&JBD_GeffSpline, loga, Geff, gg.npts);
  Create_GSL_Spline(&JBD_HSpline,    loga, H,    gg.npts);
  Create_GSL_Spline(&JBD_dHdaSpline, loga, dHda, gg.npts);

  // Free up memory
  my_free(loga);
  my_free(Geff);
  my_free(H);
  my_free(dHda);
  my_free(gg.x_arr);
  my_free(gg.y_arr);
  my_free(gg.dy_arr);
}

