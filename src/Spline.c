//  Spline routines written by Hans Winther (ICG Portsmouth) March 2017
#include "Spline.h"
#define MIN(x,y) ((x)>(y) ? (y) : (x))

// Print some info
void print_spline_info(Spline *s){
  int n = s->n;
  printf("n = %i x_start = %f  x_end = %f  type = %i\n", s->n, s->x_start, s->x_end, s->type);
  printf("x  = %f %f ... %f %f\n", s->x[0],  s->x[1],  s->x[n-2],  s->x[n-1]);
  printf("y  = %f %f ... %f %f\n", s->y[0],  s->y[1],  s->y[n-2],  s->y[n-1]);
  printf("y2 = %f %f ... %f %f\n", s->y2[0], s->y2[1], s->y2[n-2], s->y2[n-1]);
}

// Make a spline
void create_spline(Spline *s, double *x, double *y, int n, double yp1, double ypn, int type){
  double sig, p, *u, un;
  int i;

  // Allocate memory
  u     = (double *) malloc(n * sizeof(double));
  s->y2 = (double *) malloc(n * sizeof(double));
  s->y  = (double *) malloc(n * sizeof(double));
  s->x  = (double *) malloc(n * sizeof(double));
  s->arrays_allocated = 1;
  for (i = 0; i<n; i++) {
    s->y[i]  = y[i];
    s->x[i]  = x[i];
    s->y2[i] = 0.0;
  }
  
  double *y2 = s->y2;
  s->type = type;
  s->n = n;
  s->x_start = x[0];
  s->x_end = x[n-1];

  // Boundary conditions for the spline at left end
  if (yp1 > 0.99 * BC_NATURAL_SPLINE){
    y2[0] = u[0]  = 0.0;
  } else {
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  // Create spline by solving recurence relation
  for (i=1;i<n-1;i++){
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  // Boundary condition for the spline at right end
  if (ypn > 0.99 * BC_NATURAL_SPLINE){
    y2[n-1] = 0.0;
  } else {
    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    y2[n-1] = (un-0.5*u[n-2])/(0.5*y2[n-2]+1.0);
  }

  // Calculate y''
  for (i = n-2; i>=0; i--) y2[i] = y2[i]*y2[i+1]+u[i];

  // Free up memory
  free(u);
}

// Get the value from a spline
double spline_lookup(Spline *s, double x0){
  int klo, khi, k;
  double h, b, a, result;

  double *x = s->x;
  double *y = s->y;
  double *y2 = s->y2;
  double x_end = s->x_end;
  double x_start = s->x_start;
  int type = s->type;
  int n = s->n;

  if(x0 > x_end){
    printf("Warning spline_lookup: x0 = %f > x_end = %f. Returning y[x_end] = %f\n", x0, x_end, y[n-1]);
    print_spline_info(s);
    return y[n-1];
  }
  if(x0 < x_start){
    printf("Warning spline_lookup: x0 = %f < x_start = %f. Returning y[x_start] = %f\n", x0, x_start, y[0]);
    print_spline_info(s);
    return y[0];
  }

  // Find closest x-value to x0
  if (type == 1) {
    klo = MIN( (int)((x0 - x_start)/(x_end-x_start)*(n-1)) ,n-2);
    khi = klo + 1;
  } else if (type == 2) {
    klo = MIN( (int)((log(x0/x_start))/log(x_end/x_start)*(n-1)),n-2);
    khi = klo + 1;
  } else {
    klo = 0;
    khi = n-1;
    while (khi-klo>1) {
      k = (khi+klo) >> 1;
      if (x[k]>x0) {
        khi = k;
      }else {
        klo = k;
      }
    }
  }
  h = x[khi]-x[klo];
  if (h == 0.0 || khi >= n){
    printf("Spline error: h = 0 or khi >= n. The x-values muste be distinct [%i vs %i]\n", khi, n);
    print_spline_info(s);
    exit(1);
  }

  a = (x[khi]-x0)/h;
  b = (x0-x[klo])/h;
  result = (a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0);
  return result;
}

// Get the derivative of the quantity splines
double spline_lookup_dy(Spline *s, double x0){
  int klo, khi, k;
  double h, b, a, result;
  
  double *x = s->x;
  double *y = s->y;
  double *y2 = s->y2;
  double x_end = s->x_end;
  double x_start = s->x_start;
  int type = s->type;
  int n = s->n;

  if(x0 > x_end){
    printf("Warning spline_lookup_dy: x0 = %f > x_end = %f. Returning y[x_end] = %f\n", x0, x_end, y[n-1]);
    print_spline_info(s);
    return y[n-1];
  }
  if(x0 < x_start){
    printf("Warning spline_lookup_dy: x0 = %f < x_start = %f. Returning y[x_start] = %f\n", x0, x_start, y[0]);
    print_spline_info(s);
    return y[0];
  }

  // Find closest x-value to x0
  if (type == 1) {
    klo = MIN( (int)((x0 - x_start)/(x_end-x_start)*(n-1)),n-2);
    khi = klo + 1;
  } else if (type == 2) {
    klo = MIN( (int)((log(x0/x_start))/log(x_end/x_start)*(n-1)),n-2);
    khi = klo + 1;
  } else {
    // Binary search
    klo = 0;
    khi = n-1;
    while (khi-klo>1) {
      k = (khi+klo) >> 1;
      if (x[k]>x0) {
        khi = k;
      }else {
        klo = k;
      }
    }
  }
  h = x[khi]-x[klo];
  if (h == 0.0){
    printf("Spline error: h = 0. The x-values must be distinct. Exiting\n");
    print_spline_info(s);
    exit(1);
  }
  a = (x[khi]-x0)/h;
  b = (x0-x[klo])/h;
  result = (y[khi]-y[klo])/h + h/6.0*(-(3*a*a-1)*y2[klo] + (3*b*b-1)*y2[khi]);
  return result;
}

// Delete all arrays in the spline
void free_spline(Spline *s){
  if(s->arrays_allocated == 1){
    free(s->x);
    free(s->y);
    free(s->y2);
  }
  s->arrays_allocated = 0;
}

//==============================================================
// Linear interpolation in a spline-array
//==============================================================
double interpolate_from_splinearray(SplineArray *f, double k, double x){
  // First compute nearby-k values
  double kmin = f->kmin;
  double kmax = f->kmax;
  int nk      = f->nk;

  // Index to splines between current k-value
  int indlow = (int)(log(k/kmin) / log(kmax/kmin) * (nk-1));
  int indhigh = indlow + 1;

  // Bounds-check
  if(k <= kmin) {
    return spline_lookup(f->splinearray[0], x);
  }
  if(k >= kmax) {
    return spline_lookup(f->splinearray[nk-1], x);
  }

  // Linear interpolation. If out of bounds use closest value!
  if(indlow < 0.0){
    return spline_lookup(f->splinearray[0], x);
  } else if(indhigh >= nk){
    return spline_lookup(f->splinearray[nk-1], x);
  } else {
    double w = ( log(k/kmin) - log(kmax/kmin) * indlow / (nk-1) ) / ( log(kmax/kmin) / (double)(nk-1) );

    // Check that linear weight is sensible
    if(w > 1.001 || w < -0.001){
      printf("Error in interpolate spline il = %i  ih = %i n = %i  w = %f\n", indlow, indhigh, nk, w);
      printf("kmin %f  kmax %f  know %f\n", kmin, kmax, k);
      exit(1);
    }
    return spline_lookup(f->splinearray[indlow], x) * (1.0 - w) +  spline_lookup(f->splinearray[indhigh], x) * w;
  }
}

// Free up memory from a spline-container
void free_SplineArray(SplineArray *f){
  if(f->is_created){
    for(int i = 0; i < f->nk; i++){
      free_spline(f->splinearray[i]);
      free(f->splinearray[i]);
    }
    free(f->splinearray); 
  }
}

// Interpolate on a 3D grid
// Right now we don't use Splines here, but this can be implemented
// This routine assumes that we have linear spacing in the three directions
// pos is [k,k1,cosphi]
double interpolate_on_InterpolationGrid(InterpolationGrid *g, double *pos, double a, int type){
  int ix[3], ixp[3];
  double dx[3];
  int *n = g->n;

  // Compute the current time
  int ntime = (int)((log(a) - g->logamin) / (g->logamax - g->logamin) * g->ntime);
  if(ntime >= g->ntime) ntime = g->ntime-1;

  // Compute grid indices
  for(int i = 0; i < 3; i++){
    ix[i]  = (int) ( (pos[i] - g->xlow[i]) / (g->xhigh[i] - g->xlow[i]) * g->n[i] );
    ixp[i] = ix[i]+1;
    double xx = g->xlow[i] + (g->xhigh[i] - g->xlow[i]) * ix[i]/(double)(n[i]-1);
    dx[i]  = (pos[i] - xx) / (g->xhigh[i] - g->xlow[i]);
  }
  
  // Pointer to relevant grid
  double *y = NULL;
  if(type == 0) y = g->D2[ntime];
  if(type == 1) y = g->dD2dy[ntime];
  if(type == 2) y = g->ddD2ddy[ntime];
  
  // Bounds check
  for(int i = 0; i < 3; i++){
    if(ix[i] >= n[i]-1){
      ix[i]  = n[i]-1;
      ixp[i] = n[i]-1;
      dx[i]  = 0.0;
    }
    if(ix[i] < 0){
      ix[i]  = 0;
      ixp[i] = 0;
      dx[i]  = 0.0;
    }
  }

  double result;
#define TRILINEAR
#ifndef TRILINEAR

  // Fetch closest value
  int index = ix[0] + n[0]*(ix[1] + n[1]*ix[2]);
  result = y[index];

#else
  
  // Trilinear interpolation on the grid
  int ind000 = ix[0]  + n[0]*(ix[1]  + n[1]*ix[2]);
  int ind001 = ix[0]  + n[0]*(ix[1]  + n[1]*ixp[2]);
  int ind010 = ix[0]  + n[0]*(ixp[1] + n[1]*ix[2]);
  int ind100 = ixp[0] + n[0]*(ix[1]  + n[1]*ix[2]);
  int ind011 = ix[0]  + n[0]*(ixp[1] + n[1]*ixp[2]);
  int ind101 = ixp[0] + n[0]*(ix[1]  + n[1]*ixp[2]);
  int ind110 = ixp[0] + n[0]*(ixp[1] + n[1]*ix[2]);
  int ind111 = ixp[0] + n[0]*(ixp[1] + n[1]*ixp[2]);
  
  double c00 = y[ind000]*(1-dx[0]) + y[ind100] * dx[0];
  double c01 = y[ind001]*(1-dx[0]) + y[ind101] * dx[0];
  double c10 = y[ind010]*(1-dx[0]) + y[ind110] * dx[0];
  double c11 = y[ind011]*(1-dx[0]) + y[ind111] * dx[0];

  double c0 = c00 * (1-dx[1]) + c10 * dx[1];
  double c1 = c01 * (1-dx[1]) + c11 * dx[1];

  result = c0*(1-dx[2]) + c1*dx[2];

#endif

  return result;
}

