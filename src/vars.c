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
MPI_Datatype PartDataMPIType;

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

#ifdef COMPUTE_POFK

int    pofk_compute_every_step;  // Flag to turn on this option [1]: compute every step and output; [0]: don't compute at all
int    pofk_nbins;               // Number of bins in P(k) evaluation
int    pofk_bintype;             // Which binning [0]: linear [1]: log spacing 
int    pofk_subtract_shotnoise;  // Subtract shotnoise
double pofk_kmin;                // The minimum k-value in h/Mpc (should be >= 2pi/Box)
double pofk_kmax;                // The maximum k-value in h/Mpc (should be <= 2pi/Box * Nmesh)

int    pofk_compute_rsd_pofk;    // Flag to turn on computing P0,P2,P4 RSD multipole spectra. Average over 2 axes.
                                 // Compute every step: [1], Compute when outputting [2], Don't compute [0]
#endif

//===================================================
// Modified gravity variables
//===================================================

int modified_gravity_active;      // Main flag [1] is MG [0] is LCDM
int include_screening;            // 1 to include screening, 0 to keep it linear
double aexp_global;               // Global copy of current value of scale factor
int use_lcdm_growth_factors;      // If this is active the growth-factors are those in LCDM
int input_pofk_is_for_lcdm;       // If this is active then the input P(k) is rescaled according
                                  // to MG growth-factors. Useful for MG models which are LCDM at early times
int input_sigma8_is_for_lcdm;     // If this is active the sigma8 in the parameter-file is assumed to be for LCDM
                                  // Only active if input_pofk_is_for_lcdm is set. Useful to make runs with the same IC
int allocate_mg_arrays;           // This is usually 1 if we run with MG, however for some cases we don't need the
                                  // extra arrays, e.g. for simple Geff(a) models, so it's useful to have a flag for this

#ifdef MASSIVE_NEUTRINOS
int    nu_include_massive_neutrinos;
double nu_SumMassNuEV;
char   nu_FilenameTransferInfofile[500];
double OmegaNu;
double OmegaCDM;
#endif

#if defined(EQUATIONOFSTATE_PARAMETRIZATION)
double w_0; // Background parametrization
double w_a; // w(a) = w0 + wa(1-a)
#endif

#if defined(FOFRGRAVITY)
double fofr0;                    // Hu-Sawicky f(R) parameters: f(R0)            
double nfofr;                    // Hu-Sawicky f(R) parameters: n                
#elif defined(MBETAMODEL)
double beta_symm;                // Coupling beta(a->infty)
double assb_symm;                // Symmetry-breaking redshift
double range_symm;               // Range of force in Mpc/h
#elif defined(DGPGRAVITY)
double Rsmooth_global;           // Smoothing radius for density field (DGP relevant)
double rcH0_DGP;                 // DGP cross-over scale in units of c/H0
#elif defined(BRANSDICKE)
double wBD;                      // The Brans-Dicke parameter
double Omegah2;                  // The physical matter density parameter Omega * h^2
double Omegarh2;                 // The physical radiation density parameter Omegar * h^2
double Omegavh2;                 // The physical dark energy density parameter Omegav * h^2
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

// Contains the initial density and 2LPT equivalent in fourier space
// Availiable after IC generation till program exit
complex_kind *cdelta_cdm;
complex_kind *cdelta_cdm2;

// Temporary grids needed to compute growth-factors
// Allocated and deallocated as we go along
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

inline int mymod(int i, int N){
  int res = i % N;
  if(res < 0) res += N;
  return res;
}

#ifdef MATCHMAKER_HALOFINDER
int mm_run_matchmaker;     // Run matchmaker every time we output 
int mm_output_pernode;     // One file (0) or one per task (1)
int mm_output_format;      // Ascii (0), binary (2), fits (1). Fits require cfitsio library + define MATCHMAKER_USEFITS
int mm_min_npart_halo;     // Minimum number of particles per halo (20)
double mm_alloc_factor;    // Alloc factor. How much extra memory per node (1.25).
double mm_linking_length;  // FoF linking length (0.2)
double mm_dx_extra_mpc;    // Buffer radius for FoF search in units of the boxsize (3.0 Mpc/h)
#endif

// Global book-keeping variables for timesteps
int NoutputStart_global;
int timeStep_global;

