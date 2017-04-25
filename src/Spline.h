
//====================================================
// These routines are no longer used in the code...
//====================================================

#ifndef SPLINEHEADER_INC
#define SPLINEHEADER_INC
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This file contains methods for making a qubic spline of an array of (x,y)-values
// Spline routines written by Hans Winther (ICG Portsmouth) March 2017

#define ARBITRARY_SPACED_SPLINE 0
#define LINEAR_SPACED_SPLINE    1
#define LOG_SPACED_SPLINE       2
#define BC_NATURAL_SPLINE       1e30

struct SplineContainer {
    
    double *y;              // y-values
    double *x;              // x-values
    double *y2;             // y''-values computed by the spline methods
    double x_start;         // The start value in the x-array
    double x_end;           // The end value in the x-array (assuming > x_start)
    int n;                  // Number of points
    int type;               // Type of x-array (for faster lookup)
                            // NB: Wrong type will give wrong results! Use type = 0 if unsure as this will always
                            // give correct results
                            // type = 0 : arbritrary x array
                            // type = 1 : x is linear
                            // type = 2 : x is logaritmic
    int arrays_allocated;   // Book-keeping: are the arrays above allocated?
};

typedef struct SplineContainer Spline;

// Create a spline
void create_spline(Spline *s, double *x, double *y, int n, double yp1, double ypn, int type);

// Return y(x0)
double spline_lookup(Spline *s, double x0);

// Return y'(x0)
double spline_lookup_dy(Spline *s, double x0);

// Free memory
void free_spline(Spline *s);

// Print some info
void print_spline_info(Spline *s);

struct SplineArrayContainer {

  //========================================================================
  // For splines of a 2D function f(k, x) 
  // The k-values are assumed located at
  // k_i = exp[ log(kmin) + log(kmin/kmax) * i / (nk-1)]
  // This should be generalized...
  //
  // Call interpolate_from_splinearray(SplineArray *f, double k, double x)
  // to get f(k, x) after it has been computed
  //========================================================================

  double kmin;
  double kmax;
  int nk;
  int is_created;
  Spline ** splinearray;
};
typedef struct SplineArrayContainer SplineArray;

// Interpolate in a spline-array
double interpolate_from_splinearray(SplineArray *f, double k, double x);

// Free an allocated spline array
void free_SplineArray(SplineArray *f);

// For trilinear interpolation
struct InterpolationGridContainer{
  int n[3];
  int ntot;

  double logamin, logamax;
  int ntime;

  double xlow[3];
  double xhigh[3];

  int is_created;

  double **D2;
  double **dD2dy;
  double **ddD2ddy;
};
typedef struct InterpolationGridContainer InterpolationGrid;
  
// Fetch values
double interpolate_on_InterpolationGrid(InterpolationGrid *g, double *pos, double a, int type);

#endif
