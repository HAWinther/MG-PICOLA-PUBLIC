#include "wrappers.h"
#define MAX_ALLOCATIONS 1000000

int num_splines_allocated = 0;

//=================================
// Memory monitoring
//=================================
struct mapdata *memmap;
int memory_monitoring_active   =  1; 
int total_allocs               = -1;

long long total_bytes_allocated          = 0;
long long total_bytes_allocated_alltasks = 0;
long long max_bytes_allocated            = 0;
long long max_bytes_allocated_alltasks   = 0;
int times_allocated                      = 0;
int times_allocated_alltasks             = 0;
int times_freed                          = 0;
int times_freed_alltasks                 = 0;

//==============================================================================
// FFTW Wrappers to avoid having SINGLE_PRECISION-defines messing up the code
//==============================================================================

#ifdef SINGLE_PRECISION

inline plan_kind my_fftw_mpi_plan_dft_r2c_3d(int nx, int ny, int nz, 
    float_kind *regrid, complex_kind *imgrid, MPI_Comm comm, unsigned flags){
  return fftwf_mpi_plan_dft_r2c_3d((ptrdiff_t)nx, (ptrdiff_t)ny, (ptrdiff_t)nz, regrid, imgrid, comm, flags);
}

inline plan_kind my_fftw_mpi_plan_dft_c2r_3d(int nx, int ny, int nz, 
    complex_kind *imgrid, float_kind *regrid, MPI_Comm comm, unsigned flags){
  return fftwf_mpi_plan_dft_c2r_3d((ptrdiff_t)nx, (ptrdiff_t)ny, (ptrdiff_t)nz, imgrid, regrid, comm, flags);
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

inline ptrdiff_t my_fftw_mpi_local_size_3d(int nx, int ny, int nz, 
    MPI_Comm comm, ptrdiff_t *locnx, ptrdiff_t *locxstart){
  return fftwf_mpi_local_size_3d((ptrdiff_t)nx, (ptrdiff_t)ny, (ptrdiff_t)nz, comm, locnx, locxstart);
}

#else

inline plan_kind my_fftw_mpi_plan_dft_r2c_3d(int nx, int ny, int nz, 
    float_kind *regrid, complex_kind *imgrid, MPI_Comm comm, unsigned flags){
  return fftw_mpi_plan_dft_r2c_3d((ptrdiff_t)nx, (ptrdiff_t)ny, (ptrdiff_t)nz, regrid, imgrid, comm, flags);
}

inline plan_kind my_fftw_mpi_plan_dft_c2r_3d(int nx, int ny, int nz, 
    complex_kind *imgrid, float_kind *regrid, MPI_Comm comm, unsigned flags){
  return fftw_mpi_plan_dft_c2r_3d((ptrdiff_t)nx, (ptrdiff_t)ny, (ptrdiff_t)nz, imgrid, regrid, comm, flags);
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

inline ptrdiff_t my_fftw_mpi_local_size_3d(int nx, int ny, int nz, 
    MPI_Comm comm, ptrdiff_t *locnx, ptrdiff_t *locxstart){
  return fftw_mpi_local_size_3d((ptrdiff_t)nx, (ptrdiff_t)ny, (ptrdiff_t)nz, comm, locnx, locxstart);
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
  num_splines_allocated++;
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
    num_splines_allocated--;
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
  num_splines_allocated++;
}

void Set_GSL_Spline_Array(double val, double *y, int ix){
  y[ix] = val;
}

void Free_GSL_Spline(GSL_Spline *splinecontainer){
  if(splinecontainer->allocated){
    gsl_interp_accel_free(splinecontainer->xacc);
    gsl_spline_free(splinecontainer->spline);
    splinecontainer->allocated = 0;
    num_splines_allocated--;
  }
}

double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x){
  double xx = x;
  if(x < splinecontainer->xmin) xx = splinecontainer->xmin; 
  if(x > splinecontainer->xmax) xx = splinecontainer->xmax; 
  return gsl_spline_eval(splinecontainer->spline, xx, splinecontainer->xacc);
}

//=================================
// Wrapper on malloc
//=================================
inline void *my_malloc(size_t bytes){
  void *ptr =  malloc(bytes);
  if(memory_monitoring_active){
    size_t addr = (size_t) ptr;
    add_to_memmap(addr, bytes);
    total_bytes_allocated += (long long) bytes;
    if(total_bytes_allocated > max_bytes_allocated){
      max_bytes_allocated = total_bytes_allocated;
    }
  }
  times_allocated++;
  return ptr;
}

//=================================
// Wrapper on calloc
//=================================
inline void* my_calloc(size_t n, size_t bytes_per_element){
  void *ptr = calloc(n, bytes_per_element);
  size_t bytes = n * bytes_per_element;
  if(memory_monitoring_active){
    size_t addr = (size_t) ptr;
    add_to_memmap(addr, bytes);
    total_bytes_allocated += (long long) bytes;
    if(total_bytes_allocated > max_bytes_allocated){
      max_bytes_allocated = total_bytes_allocated;
    }
  }
  times_allocated++;
  return ptr;
} 

//======================================================
// Wrapper on free
// We should add a proper hashmap or similar to store the
// (addr,bytes) pairs below
//======================================================
inline void my_free(void *ptr){
  if(memory_monitoring_active){
    size_t addr = (size_t)ptr;
    int i = find_index_in_memmap(addr);
    if(i < 0){
      printf("Warning: my_free is freeing something not allocated by my_malloc\n");
    } else {
      size_t bytes = memmap[i].bytes;
      total_bytes_allocated -= bytes;
      memmap[i].addr  = 0;
      memmap[i].bytes = 0;
    }
  }
  times_freed++;
  free(ptr);
}

void init_memory_monitoring(){
  if(memory_monitoring_active){
    memmap = malloc(sizeof(struct mapdata) * MAX_ALLOCATIONS);
    total_allocs = 0;
    total_bytes_allocated = 0;
    max_bytes_allocated = 0;
  }
}

void cleanup_memory_monitor(){
  if(memory_monitoring_active){
    free(memmap);
    total_allocs = -1;
  }
}

int find_index_in_memmap(size_t addr){
  for(int i = 0; i < total_allocs; i++){
    if(memmap[i].addr == addr) return i;
  }
  return -1;
}

void add_to_memmap(size_t addr, size_t bytes){
  if(total_allocs == -1) init_memory_monitoring();
  if(total_allocs >= MAX_ALLOCATIONS){
    printf("Warning: total allocs exceeds max allocations. Memory monitoring will fail to work\n");
    free(memmap);
    memory_monitoring_active = 0;
    return;
  }
  memmap[total_allocs].addr = addr;
  memmap[total_allocs].bytes = bytes;
  total_allocs++;
}

void print_memory_thistask(){
  printf("\n===============================\n"); 
  printf(  "Memory summary Task[%i]:       \n", ThisTask);
  printf(  "===============================\n"); 
  printf("Memory allocated now: [%5.3f] MB (max over tasks)\n", (double) total_bytes_allocated / (1024. * 1024.));
  printf("Peak memory usage   : [%5.3f] MB (max over tasks)\n", (double) max_bytes_allocated / (1024. * 1024.));
  printf("Allocated %i times and freed %i times (sum over tasks) nSplines [%i]\n", times_allocated, times_freed, num_splines_allocated);
  printf("\n");
  fflush(stdout);
}

void print_memory_summary(){
  
  MPI_Allreduce(&total_bytes_allocated, &total_bytes_allocated_alltasks, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_bytes_allocated,   &max_bytes_allocated_alltasks,   1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&times_allocated,       &times_allocated_alltasks,       1, MPI_INT,       MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&times_freed,           &times_freed_alltasks,           1, MPI_INT,       MPI_SUM, MPI_COMM_WORLD);
    
  if(ThisTask == 0){
    printf("\n===============================\n"); 
    printf(  "Memory summary:                \n");
    printf(  "===============================\n"); 
    printf("Memory allocated now: [%5.3f] MB (max over tasks)\n", (double) total_bytes_allocated_alltasks / (1024. * 1024.));
    printf("Peak memory usage   : [%5.3f] MB (max over tasks)\n", (double) max_bytes_allocated_alltasks / (1024. * 1024.));
    printf("We have allocated %i times and freed %i times (sum over tasks)\n", times_allocated_alltasks, times_freed_alltasks);
    printf("Splines allocated: %i\n",num_splines_allocated);
    printf("\n");
    fflush(stdout);
  }
}

