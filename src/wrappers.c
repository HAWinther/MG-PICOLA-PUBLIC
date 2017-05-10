#include "wrappers.h"

//==============================================================================
// FFTW Wrappers to avoid having SINGLE_PRECISION-defines messing up the code
//==============================================================================

#ifdef SINGLE_PRECISION

inline plan_kind my_fftw_mpi_plan_dft_r2c_3d(int nx, int ny, int nz, float_kind *regrid, complex_kind *imgrid, MPI_Comm comm, unsigned flags){
  return fftwf_mpi_plan_dft_r2c_3d(nx, ny, nz, regrid, imgrid, comm, flags);
}

inline plan_kind my_fftw_mpi_plan_dft_c2r_3d(int nx, int ny, int nz, complex_kind *imgrid, float_kind *regrid, MPI_Comm comm, unsigned flags){
  return fftwf_mpi_plan_dft_c2r_3d(nx, ny, nz, imgrid, regrid, comm, flags);
}

inline void my_fftw_destroy_plan(plan_kind fftwplan){
  fftwf_destroy_plan(fftwplan);
}

inline void my_fftw_execute(plan_kind fftwplan){
  timer_start(_FFT);
  fftwf_execute(fftwplan);
  timer_stop(_FFT);
}

inline void my_fftw_mpi_cleanup(){
  fftwf_mpi_cleanup();
}

inline void my_fftw_mpi_init(){
  fftwf_mpi_init();
}

inline ptrdiff_t my_fftw_mpi_local_size_3d(int nx, int ny, int nz, MPI_Comm comm, ptrdiff_t *locnx, ptrdiff_t *locxstart){
  return fftwf_mpi_local_size_3d(nx, ny, nz, comm, locnx, locxstart);
}

#else

inline plan_kind my_fftw_mpi_plan_dft_r2c_3d(int nx, int ny, int nz, float_kind *regrid, complex_kind *imgrid, MPI_Comm comm, unsigned flags){
  return fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, regrid, imgrid, comm, flags);
}

inline plan_kind my_fftw_mpi_plan_dft_c2r_3d(int nx, int ny, int nz, complex_kind *imgrid, float_kind *regrid, MPI_Comm comm, unsigned flags){
  return fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, imgrid, regrid, comm, flags);
}

inline void my_fftw_destroy_plan(fftw_plan fftwplan){
  fftw_destroy_plan(fftwplan);
}

inline void my_fftw_execute(fftw_plan fftwplan){
  timer_start(_FFT);
  fftw_execute(fftwplan);
  timer_stop(_FFT);
}

inline void my_fftw_mpi_cleanup(){
  fftw_mpi_cleanup();
}

inline void my_fftw_mpi_init(){
  fftw_mpi_init();
}

inline ptrdiff_t my_fftw_mpi_local_size_3d(int nx, int ny, int nz, MPI_Comm comm, ptrdiff_t *locnx, ptrdiff_t *locxstart){
  return fftw_mpi_local_size_3d(nx, ny, nz, comm, locnx, locxstart);
}

#endif

//==============================================================================
// GSL wrappers
//==============================================================================

void Create_GSL_2D_Spline(GSL_2D_Spline *splinecontainer, double *x, double *y, double *z, int nx, int ny){
  splinecontainer->T      = gsl_interp2d_bicubic;
  splinecontainer->spline = gsl_spline2d_alloc(splinecontainer->T, nx, ny);
  splinecontainer->xacc   = gsl_interp_accel_alloc();
  splinecontainer->yacc   = gsl_interp_accel_alloc();
  splinecontainer->xmin   = x[0];
  splinecontainer->xmax   = x[nx-1];
  splinecontainer->ymin   = y[0];
  splinecontainer->ymax   = y[ny-1];
  gsl_spline2d_init(splinecontainer->spline, x, y, z, nx, ny);
  splinecontainer->allocated = 1;
}

void Set_GSL_2D_Spline_Array(double val, double *z, int nx, int ix, int iy){
  z[iy * nx + ix] = val;
}

void Free_GSL_2D_Spline(GSL_2D_Spline *splinecontainer){
  if(splinecontainer->allocated){
    gsl_interp_accel_free(splinecontainer->xacc);
    gsl_interp_accel_free(splinecontainer->yacc);
    gsl_spline2d_free(splinecontainer->spline);
    splinecontainer->allocated = 0;
  }
}

double Lookup_GSL_2D_Spline(GSL_2D_Spline *splinecontainer, double x, double y){
  double xx = x, yy = y;
  if(x < splinecontainer->xmin) xx = splinecontainer->xmin; 
  if(x > splinecontainer->xmax) xx = splinecontainer->xmax; 
  if(y < splinecontainer->ymin) yy = splinecontainer->ymin; 
  if(y > splinecontainer->ymax) yy = splinecontainer->ymax; 
  return gsl_spline2d_eval(splinecontainer->spline, xx, yy, splinecontainer->xacc, splinecontainer->yacc);
}

void Create_GSL_Spline(GSL_Spline *splinecontainer, double *x, double *y, int nx){
  splinecontainer->xacc   = gsl_interp_accel_alloc();
  splinecontainer->xmin   = x[0];
  splinecontainer->xmax   = x[nx-1];
  splinecontainer->spline = gsl_spline_alloc(gsl_interp_cspline, nx);
  gsl_spline_init(splinecontainer->spline, x, y, nx);
  splinecontainer->allocated = 1;
}

void Set_GSL_Spline_Array(double val, double *y, int ix){
  y[ix] = val;
}

void Free_GSL_Spline(GSL_Spline *splinecontainer){
  if(splinecontainer->allocated){
    gsl_interp_accel_free(splinecontainer->xacc);
    gsl_spline_free(splinecontainer->spline);
    splinecontainer->allocated = 0;
  }
}

double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x){
  double xx = x;
  if(x < splinecontainer->xmin) xx = splinecontainer->xmin; 
  if(x > splinecontainer->xmax) xx = splinecontainer->xmax; 
  return gsl_spline_eval(splinecontainer->spline, xx, splinecontainer->xacc);
}

