#ifndef INCWRAPPERS
#define INCWRAPPERS
#include "timer.h"
#include "vars.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

//========================================
//
// This file contains several wrappers
// for FFTW and GSL to avoid cluttering
// up the code with defines etc. and
// to ease the creation of new objects
//
// We also have wrappers for malloc
// and free to keep track of how much
// memory we have allocated. If allocating
// with my_alloc one must free with my_free
// for this to work.
//
//========================================

plan_kind my_fftw_mpi_plan_dft_r2c_3d(int nx, int ny, int nz, 
    float_kind *regrid, complex_kind *imgrid, MPI_Comm comm, unsigned flags);
plan_kind my_fftw_mpi_plan_dft_c2r_3d(int nx, int ny, int nz, 
    complex_kind *imgrid, float_kind *regrid, MPI_Comm comm, unsigned flags);
void      my_fftw_destroy_plan(plan_kind fftwplan);
void      my_fftw_execute(plan_kind fftwplan);
void      my_fftw_mpi_cleanup();
void      my_fftw_mpi_init();
ptrdiff_t my_fftw_mpi_local_size_3d(int nx, int ny, int nz, 
    MPI_Comm comm, ptrdiff_t *locnx, ptrdiff_t *locxstart);

typedef struct GSL_Spline {
  gsl_spline *spline;
  gsl_interp_accel *xacc;
  double xmin, xmax;
  int allocated;
} GSL_Spline;

typedef struct GSL_2D_Spline {
  const gsl_interp2d_type *T;
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  double xmin, xmax;
  double ymin, ymax;
  int allocated;
} GSL_2D_Spline;

void   Create_GSL_Spline(GSL_Spline *splinecontainer, double *x, double *y, int nx);
void   Free_GSL_Spline(GSL_Spline *splinecontainer);
double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x);
void   Set_GSL_Spline_Array(double val, double *y, int ix);

void   Create_GSL_2D_Spline(GSL_2D_Spline *splinecontainer, double *x, double *y, double *z, int nx, int ny);
void   Free_GSL_2D_Spline(GSL_2D_Spline *splinecontainer);
double Lookup_GSL_2D_Spline(GSL_2D_Spline *splinecontainer, double x, double y);
void   Set_GSL_2D_Spline_Array(double val, double *z, int nx, int ix, int iy);

struct mapdata {
  size_t addr;
  size_t bytes;
};

void* my_calloc(size_t n, size_t bytes_per_element);
void* my_malloc(size_t bytes);
void  my_free(void *ptr);
void init_memory_monitoring();
void cleanup_memory_monitor();
int find_index_in_memmap(size_t addr);
void add_to_memmap(size_t addr, size_t bytes);
void print_memory_summary();

#endif
