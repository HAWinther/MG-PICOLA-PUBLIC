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

//======================================================================================//
// This file contains the initialistion of all the external variables in the header file//
//======================================================================================//

#include "vars.h"
#include "msg.h"
#include "timer.h"

//===================================================
// MPI variables
//===================================================
int ierr;             // The return value for mpi routines
int NTask;            // The total number of tasks
int ThisTask;         // The rank of each task
int LeftTask;         // The first neighbouring task on the left containing particles
int RightTask;        // The first neighbouring task on the right containing particles
MPI_Status status;    // The MPI error status
MPI_Request request;  // The continue directive for non-blocking sends

//===================================================
// Global variables for the grids
//===================================================
int NTaskWithN;           // The number of tasks that actually have particles
int last_slice;           // The last slice of the density/force grids (maybe equal to alloc_local)
int * Slab_to_task;       // The task to which each slice is assigned
int * Part_to_task;       // The task to which each particle position is assigned
int * Local_nx_table;     // The number of slices on each of the tasks
int * Local_np_table;     // The number of particle grid slices on each of the tasks
float_kind * N11;         // The force grid in the X direction
float_kind * N12;         // The force grid in the Y direction
float_kind * N13;         // The force grid in the Z direction
float_kind * density;     // The density grid
ptrdiff_t Local_nx;       // The number of slices on the task
ptrdiff_t Local_np;       // The number of particle grid slices on the task
ptrdiff_t Total_size;     // The total byte-size of the grids on each processor
ptrdiff_t alloc_local;    // The byte-size returned by FFTW required to allocate the density/force grids
ptrdiff_t alloc_slice;    // The byte-size of a slice of the density/force grids
ptrdiff_t Local_x_start;  // The global start of the slices on the task
ptrdiff_t Local_p_start;  // The global start of the particle grid slices on the task
complex_kind * P3D;       // Pointer to the complex, FFT'ed density grid (use in-place FFT)
complex_kind * FN11;      // Pointer to the complex, FFT'ed N11 force grid (use in-place FFT)
complex_kind * FN12;      // Pointer to the complex, FFT'ed N12 force grid (use in-place FFT)
complex_kind * FN13;      // Pointer to the complex, FFT'ed N13 force grid (use in-place FFT)
plan_kind plan;           // The plan for the in-place FFT of the density grid
plan_kind p11,p12,p13;    // Plans for the in-place FFT's of the forces grids 

//===================================================
// Modified gravity variables
//===================================================

int modified_gravity_active;      // Main flag [1] is MG [0] is LCDM
int include_screening;            // 1 to include screening, 0 to keep it linear
double aexp_global;               // Global copy of current value of scale factor
int use_lcdm_growth_factors;      // If this is active the growth-factors are those in LCDM
int input_sigma8_is_for_lcdm;     // If this is active P(k) is rescaled according to
                                  // the mg growth-factor when making IC
#if defined(FOFRGRAVITY) || defined(MBETAMODEL)
double fofr0;                     // Hu-Sawicky f(R) parameters: f(R0)            
double nfofr;                     // Hu-Sawicky f(R) parameters: n                
#elif defined(DGPGRAVITY)
double Rsmooth_global;            // Smoothing radius for density field (DGP relevant)
double rcH0_DGP;                  // DGP cross-over scale in units of c/H0
#endif

#if defined(MBETAMODEL)
Spline *phi_of_a_spline = NULL;
#endif

float_kind * mgarray_one;         // Modified gravity arrays
float_kind * mgarray_two;         // ...
complex_kind * P3D_mgarray_one;   // k-space arrays
complex_kind * P3D_mgarray_two;   // ...
plan_kind plan_mg_phinewton;      // FFT plans
plan_kind plan_mg_phik;           // ...

//===================================================
// Units
//===================================================
double G;                               // The unit-less Gravitational constant
double Light;                           // The unit-less speed of light
double Hubble;                          // The unit-less Hubble constant
double UnitMass_in_g;                   // The unit mass (in g/cm) used in the code, read in from run parameters
double UnitTime_in_s;                   // The unit time (in s) used for the code, calculated from unit length and velocity
double UnitLength_in_cm;                // The unit length (in cm/h) used in the code, read in from run parameters file
double UnitVelocity_in_cm_per_s;        // The unit velocity (in cm/s) used in the code, read in from run parameters file
double InputSpectrum_UnitLength_in_cm;  // The unit length (in cm/h) of the tabulated input spectrum, read in from run parameters

//===================================================
// Gadget-Style header
//===================================================
#ifdef GADGET_STYLE
struct io_header_1 header;
#endif

//===================================================
// Cosmological parameters (at z=0)
//===================================================
char OutputRedshiftFile[500];  // The list of output redshifts
int timeSteptot;               // The total number of timsteps made
double Fnl;                    // The primordial non-gaussianity parameter for local, equilateral or orthogonal
double Anorm;                  // The normalisation of the power spectrum/ transfer function
double Omega;                  // The total matter density, CDM+Baryon
double Sigma8;                 // The normalisation of the power spectrum 
double FnlTime;                // The scale factor at which fnl kicks in
double DstartFnl;              // The growth factor for the initial potential 
double ShapeGamma;             // The paramerisation of the Efstathiou power spectrum
double OmegaBaryon;            // The baryonic matter density
double HubbleParam;            // The normalised Hubble parameter, h=H/100
double Fnl_Redshift;           // The redshift at which the nongaussian f_nl potential is computed
double Init_Redshift;          // The redshift at which to begin timestepping
double PrimordialIndex;        // The spectral index, n_s
struct Outputs * OutputList;   // The output redshifts of the simulation and the number of steps between outputs

#ifdef GENERIC_FNL
//===================================================
// Kernels to include general non-gaussian models
//===================================================
int NKernelTable;                 // The length of the kernel lookup table
struct kern_table * KernelTable;  // The kernel lookup table
#endif

//===================================================
// Particle data and pointers
//===================================================

double sumxyz[3];
double sumDxyz[3];

#ifdef SCALEDEPENDENT

// Stored 1LPT displacment scalar in k-space at z = 0 
complex_kind *(cdisp_store[3]);
float_kind *(disp_store[3]);

// Stored 2LPT displacment scalar in k-space at z = 0 
complex_kind *(cdisp2_store[3]);
float_kind *(disp2_store[3]);

// Temporary fields needed to compute growth-factors
float_kind *ZA_D[3];
float_kind *disp_D[3];

#endif

#ifdef MEMORY_MODE

float * Disp[3];    // Vectors to hold the particle displacements each timestep
float * ZA[3];      // Vectors to hold the Zeldovich displacements before particle initialisation
float * LPT[3];     // Vectors to hold the 2LPT displacements before particle initialisation

#else

float_kind * Disp[3];    // Vectors to hold the particle displacements each timestep
float_kind * ZA[3];      // Vectors to hold the Zeldovich displacements before particle initialisation
float_kind * LPT[3];     // Vectors to hold the 2LPT displacements before particle initialisation

#endif

struct part_data *P;

//===================================================
// Simulation variables
//===================================================
char FileBase[500];             // The base output filename
char OutputDir[500];            // The output directory
int Nmesh;                      // The size of the displacement, density and force grids (in 1-D) 
int Nsample;                    // The number of particles (in 1-D)
int UseCOLA;                    // Whether or not to use the COLA modifications
int Noutputs;                   // The number of output times
int NumFilesWrittenInParallel;  // The maximum number of files to be written out in parallel
unsigned int NumPart;           // The number of particles on each processor
unsigned long long TotNumPart;  // The total number of particles in the simulation
double Box;                     // The edge length of the simulation
double Buffer;                  // The amount of extra memory of each processor to compensate for moving particles
#ifdef LIGHTCONE
int * writeflag;                // A flag to tell the code whether to write a new file or append onto an existing one.
int * repflag;                  // A flag to say whether we need to check inside a given replicate
int Nrep_neg_x;                 // The number of replicated boxes in the negative x direction
int Nrep_neg_y;                 // The number of replicated boxes in the negative y direction
int Nrep_neg_z;                 // The number of replicated boxes in the negative z direction
int Nrep_pos_x;                 // The number of replicated boxes in the positive x direction
int Nrep_pos_y;                 // The number of replicated boxes in the positive y direction
int Nrep_pos_z;                 // The number of replicated boxes in the positive z direction
int Nrep_neg_max[3];            // The maximum number of replicated boxes in the negative directions
int Nrep_pos_max[3];            // The maximum number of replicated boxes in the positive directions
unsigned int * Noutput;         // The number of particles that we output in each slice (original and replicated)
double Origin_x;                // The x-position of the lightcone origin
double Origin_y;                // The y-position of the lightcone origin
double Origin_z;                // The z-position of the lightcone origin
#endif

//===================================================
// 2LPT specific
//===================================================
char FileWithInputSpectrum[500];  // The file containing the input power spectrum
char FileWithInputTransfer[500];  // The file containing the input transfer function
char FileWithInputKernel[500];    // The file containing the input nongaussian kernel
int Seed;                         // The random seed to generate to realisation
int SphereMode;                   // Whether to use a sphere or a cube in k-space
int WhichSpectrum;                // Which power spectrum to use
int WhichTransfer;                // Which transfer function to use

//===================================================
// COLA specific
//===================================================
int fullT;        // The time dependence of the velocity, hardcoded (see README)
int StdDA;        // The time dependence of the displacement, hardcoded (see README)
double nLPT;      // Parameterisation of the time dependence of the velocity, hardcoded (see README)

//===================================================
// Read IC from RAMSES / Gadget / Ascii files
//===================================================
char InputParticleFileDir[200];
char InputParticleFilePrefix[200];
int  NumInputParticleFiles;
int  RamsesOutputNumber;
int  TypeInputParticleFiles;
int  ReadParticlesFromFile;

//==============================================================================
// FFTW Wrappers to avoid having SINGLE_PRECISION-defines messing up the code
//==============================================================================

inline plan_kind my_fftw_mpi_plan_dft_r2c_3d(int nx, int ny, int nz, float_kind *regrid, complex_kind *imgrid, MPI_Comm comm, unsigned flags){
  timer_start(_FFT);
#ifdef SINGLE_PRECISION
  return fftwf_mpi_plan_dft_r2c_3d(nx, ny, nz, regrid, imgrid, comm, flags);
#else
  return fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, regrid, imgrid, comm, flags);
#endif
  timer_stop(_FFT);
}
inline plan_kind my_fftw_mpi_plan_dft_c2r_3d(int nx, int ny, int nz, complex_kind *imgrid, float_kind *regrid, MPI_Comm comm, unsigned flags){
  timer_start(_FFT);
#ifdef SINGLE_PRECISION
  return fftwf_mpi_plan_dft_c2r_3d(nx, ny, nz, imgrid, regrid, comm, flags);
#else
  return fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, imgrid, regrid, comm, flags);
#endif
  timer_stop(_FFT);
}
inline void my_fftw_destroy_plan(fftw_plan fftwplan){
#ifdef SINGLE_PRECISION
  fftwf_destroy_plan(fftwplan);
#else
  fftw_destroy_plan(fftwplan);
#endif
}
inline void my_fftw_execute(fftw_plan fftwplan){
  timer_start(_FFT);
#ifdef SINGLE_PRECISION
  fftwf_execute(fftwplan);
#else
  fftw_execute(fftwplan);
#endif
  timer_stop(_FFT);
}
inline void my_fftw_mpi_cleanup(){
#ifdef SINGLE_PRECISION
  fftwf_mpi_cleanup();
#else
  fftw_mpi_cleanup();
#endif
}
inline void my_fftw_mpi_init(){
#ifdef SINGLE_PRECISION
  fftwf_mpi_init();
#else
  fftw_mpi_init();
#endif
}
inline ptrdiff_t my_fftw_mpi_local_size_3d(int nx, int ny, int nz, MPI_Comm comm, ptrdiff_t *locnx, ptrdiff_t *locxstart){
#ifdef SINGLE_PRECISION
  return fftwf_mpi_local_size_3d(nx, ny, nz, comm, locnx, locxstart);
#else
  return fftw_mpi_local_size_3d(nx, ny, nz, comm, locnx, locxstart);
#endif
}

inline int mymod(int i, int N){
  int res = i % N;
  if(res < 0) res += N;
  return res;
}

