#ifndef INCPROTO
#define INCPROTO

//==========================================================================//
//  Copyright (c) 2013       Cullan Howlett & Marc Manera,                  //
//                           Institute of Cosmology and Gravitation,        //
//                           University of Portsmouth.                      //
//                                                                          //
// MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017            //
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

//==========================================================================//
// v1.0: This file contains all the prototypes for function used in the code//
//==========================================================================//

#include "wrappers.h"

// Some useful functions
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define pow2(x)  ((x)*(x))
#define pow3(x)  ((x)*(x)*(x))
#define pow4(x)  ((x)*(x)*(x)*(x))
#define pow5(x)  ((x)*(x)*(x)*(x)*(x))

//===================================================
// main.c
//===================================================

void Output_Info(double A);
void Output(double A, double AF, double AFF, double dDdy, double dD2dy, double z_from_outputlist);
void Kick(double AI, double AF, double A, double Di);
void Drift(double A, double AFF, double AF, double Di, double Di2);
void write_gadget_header(FILE *fp, double A);

//===================================================
// Modified gravity routines mg.h
//===================================================

void   ComputeFifthForce();
void   ComputeFifthForce_PotentialScreening();
void   ComputeFifthForce_DensityScreening();
void   ComputeFifthForce_GradientScreening();
void   ComputeFifthForce_TimeDepGeffModels();
void   check_real_space_grid(float_kind *grid);
double smoothing_filter(double kR);
double beta_of_a(double a);
double mass2_of_a(double a);
double coupling_function(double a);
double screening_factor_potential(double a, double phinewton);
double screening_factor_density(double a, double density);
double screening_factor_gradient(double a, double DPhi2);
void   DivideByLaplacian(complex_kind *densityk, complex_kind *phinewtonk);
void   Density_to_DPhiNewtonk(complex_kind *densityk, complex_kind *DPhi_i, int axes);
void   EffDensitykToPhiofk(complex_kind* P3D_densityeffk, complex_kind *P3D_phik);
void   ComputeFifthForce();
void   ComputeFifthForceDGP();
void   SmoothDensityField(complex_kind *densityk, complex_kind *densityk_smooth, double Rsmooth);
void   CopyDensityArray();
void   AllocateMGArrays();
void   FreeMGArrays();

#ifdef BRANSDICKE
void   JBD_Solve_Background(double w, double omegamh2, double omegavh2, double omegarh2, double *h);
double JBD_Hubble_of_a(double a);
double JBD_dHubbleda_of_a(double a);
double JBD_GeffG_of_a(double a);
#endif

//===================================================
// Updated cosmo routines cosmo.c
//===================================================

void   solve_for_growth_factors();
void   free_up_splines();
void   free_CAMB_splines();

// Growth factor splines
double growth_D(double a);
double growth_dDdy(double a);
double growth_ddDddy(double a);
double growth_D2(double a);
double growth_dD2dy(double a);
double growth_ddD2ddy(double a);

// Growth factor for LCDM splines
double growth_DLCDM(double a);
double growth_dDLCDMdy(double a);
double growth_ddDLCDMddy(double a);
double growth_D2LCDM(double a);
double growth_dD2LCDMdy(double a);
double growth_ddD2LCDMddy(double a);

// Growth factor for LCDM fitting functions
double growth_D_LCDMFit(double a);
double growth_dDdy_LCDMFit(double a);
double growth_ddDddy_LCDMFit(double a);
double growth_D2_LCDMFit(double a);
double growth_dD2dy_LCDMFit(double a);
double growth_ddD2ddy_LCDMFit(double a);

// Right hand side of the ODE system for the growth-factors
int ode_first_order_growth_D(double x, const double D[], double dDdy[], void *params);
int ode_first_order_growth_DLCDM(double x, const double D[], double dDdy[], void *params);
int ode_second_order_growth_D2(double x, const double D2[], double dD2dy[], void *params);
int ode_second_order_growth_D2LCDM(double x, const double D2[], double dD2dy[], void *params);

// The second order growth ODE needs the first order growth factor plus k-value (if scale-dependent)
// and for the second order growth kernel (if we compute it)
struct ode_second_order_growth_parameters{
  double k_value;
  double k1_value;
  double k2_value;
  double costheta_value;
};

// Cosmology functions
double hubble(double a);
double Qfactor(double a);
double dhubbleda(double a);
double beta_DGP(double a);
double GeffoverG(double a, double k);
double Factor_2LPT(double a);

// For modified gravity models defined by m(a) beta(a)
double phi_of_a(double a);
void   compute_phi_of_a();
double phi_of_a_from_spline(double a);
double integrand_phiofa(double x, void *params);

double mg_pofk_ratio(double k, double a);  // Ratio of P_MG(k,a) / P_LCDM(k,a)
double mg_sigma8_enhancement(double a);    // Linear theory enhancement of sigma8
void   read_mg_parameters(void **addr, char (*tag)[50], int *id, int *nt);
void   init_modified_version();

#ifdef SCALEDEPENDENT

#define LPT_ORDER_ONE 1
#define LPT_ORDER_TWO 2
#define FIELD_D       0
#define FIELD_dDdy    1
#define FIELD_ddDddy  2
#define FIELD_deltaD  3

// Scaledependent growth-factors
void  integrate_first_order_scale_dependent_growth_ode(double *x_arr, double *d_arr, double *q_arr, 
    double know, int npts);
void  integrate_second_order_scale_dependent_growth_ode(double *x_arr, double *d_arr, double *q_arr, 
    struct ode_second_order_growth_parameters *ode_D2_params, int npts);
void   calculate_scale_dependent_growth_factor();

double growth_D_scaledependent(double k, double a);
double growth_dDdy_scaledependent(double k, double a);
double growth_ddDddy_scaledependent(double k, double a);

double growth_D2_scaledependent(double k, double a);
double growth_dD2dy_scaledependent(double k, double a);
double growth_ddD2ddy_scaledependent(double k, double a);

// 2LPT
void from_cdisp_store_to_ZA(double A, double AF, double AFF, int firststep, int LPTorder);
void assign_displacment_field_to_particles(double A, double AF, double AFF, int firststep, int LPTorder);
#endif

//===================================================
// cosmo.c
//===================================================

double gpQ(double a);
double DERgpQ(double a);
double decayD(double a);
double DprimeQ(double a);
double Qfactor_LCDM(double a);
double AofTime(double y);
double growthD(double a);
double growthD2(double a);
double growthD2v(double a);
double CosmoTime(double af);
double growthDtemp(double a);
double growthD2temp(double a);
double SqStd(double ai,double af);
double SphiStd(double ai,double af);
double fun(double x, void * params);
double funSqStd(double a, void * params);
double AofTimeFun(double a,void * params);
double Sq(double ai,double af,double aRef);
double funSphiStd(double a, void * params);
double Sphi(double ai,double af,double aRef);
double CosmoTimeFun (double a, void * params);
double GrowthFactor(double astart, double aend);

//===================================================
// auxPM.c
//===================================================

void Forces(void);
void PtoMesh(void);
void MtoParticles(void);
void MoveParticles(void);
void GetDisplacements(void);
void FatalError(char * errmsg);
#if (MEMORY_MODE || SINGLE_PRECISION)
float periodic_wrap(float x);
#else
double periodic_wrap(double x);
#endif
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);

//===================================================
// read_param.c
//===================================================

void read_outputs(void);
void read_parameterfile(char * fname);
int sort_redshift(const void * Item1, const void * Item2);

//===================================================
// kernel.c
//===================================================

#ifdef GENERIC_FNL
void read_kernel_table(void);
#endif

//===================================================
// lightcone.c
//===================================================

#ifdef LIGHTCONE
void set_lightcone(void);
void Output_Info_Lightcone(void);
void Output_Lightcone(unsigned int * pc, unsigned int blockmaxlen, float * block);
void flag_replicates(double Rcomov_old, double Rcomov_new, double boundary);
void Drift_Lightcone(double A, double AFF, double AF, double Di, double Di2);
double nearest_dist(double px, double py, double ix, double iy, double jx, double jy, double boundary);
#endif

//===================================================
// 2LPT.c
//===================================================

void set_units(void);
void initialize_ffts(void);
void initialize_parts(void);
void displacement_fields(void);

//===================================================
// power.c
//===================================================

void print_spec(void);
void free_powertable(void);
void read_power_table(void);
void free_transfertable(void);
void read_transfer_table(void);
void initialize_powerspectrum(void);
void initialize_transferfunction(void);
int compare_logk(const void *a, const void *b);
int compare_transfer_logk(const void *a, const void *b);
double fnl(double x);
double TransferFunc(double k);
double PowerSpec(double kmag);
double TopHatSigma2(double R);
double PowerSpec_EH(double k);
double TransferFunc_EH(double k);
double PowerSpec_Tabulated(double k);
double PowerSpec_Efstathiou(double k);
double TransferFunc_Tabulated(double k);
double sigma2_int(double k, void * params);

//===================================================
// Read IC from RAMSES / Gadget / Ascii file instead 
// of computing it [readICfromfile.c]
//===================================================

#ifdef READICFROMFILE
void AssignDisplacementField(complex_kind *(cdisp[3]));
void ReadFilesMakeDisplacementField(void);
int  read_gadget_file(char *filedir, char *fileprefix, int filenum, char *buffer, int *npart_read);
int  read_ascii_file(char *filedir, char *fileprefix, int filenum, char *buffer, int *npart_read);
int  read_ramses_file(char *filedir, int outnumber, int filenum, char *buffer, int *npart_read);
int  find_maxpart_gadgetfiles(char *outputdir, char *fileprefix, int nfiles);
int  find_maxpart_asciifiles(char *outputdir, char *fileprefix, int nfiles);
int  find_maxpart_ramsesfiles(char *outputdir, int outnumber, int nfiles);
void read_gadget_header(FILE *fp);
void read_ramses_header(FILE *fp);
void read_double_vec(FILE* fp, double *buffer, int n);
void read_float_vec(FILE* fp, float *buffer, int n);
void read_int_vec(FILE* fp, int *buffer, int n);
int  read_int(FILE* fp);
int  ProcessParticlesSingleFile(char *buffer, int npart_loc);
#endif

#ifdef SCALEDEPENDENT
int    ode_second_order_growth_kernel_D2(double x, const double D2[], double dD2dy[], void *params);
void   solve_for_second_order_kernel();
double fofr_pi_factor(double k, double a);
double second_order_kernel(double k, double k1, double k2, double costheta, double a);
#endif

#ifdef MASSIVE_NEUTRINOS
void read_and_spline_total_power_spectrum(char *pofkfilename);
void read_and_spline_transfer_functions(char *transferfileinfofilename);
void read_single_transfer_function_file(char *filename, double *logk, double *transfer_function_nu, 
    double *transfer_function_cdm, double *transfer_function_baryon, double *transfer_function_total, int *npts);
double get_nu_transfer_function(double k, double a);
double get_cdm_transfer_function(double k, double a);
double get_baryon_transfer_function(double k, double a);
double get_cdm_baryon_transfer_function(double k, double a);
double get_total_transfer_function(double k, double a);
double get_total_power_spectrum(double k);
double get_cdm_baryon_power_spectrum(double k);
double get_neutrino_power_spectrum(double k);
#endif

#ifdef COMPUTE_POFK
void compute_power_spectrum(complex_kind *dens_k, double a, char *label);
void adjust_pofk_parameters(int *nbins, int *bintype, int *subtract_shotnoise, double *kmin, double *kmax);
void compute_RSD_powerspectrum(double A, int dDdy_set_in_particles);
void PtoMesh_RSD(float_kind *dens, int axis, double A);
void bin_up_RSD_power_spectrum(complex_kind *dens_k, struct RSD_Pofk_Data *pofk_data);
int pofk_bin_index(double kmag, double kmin, double kmax, int nbins, int bintype);
double k_from_index(int index, double kmin, double kmax, int nbins, int bintype);
#endif

#ifdef MATCHMAKER_HALOFINDER
void MatchMaker(struct PicolaToMatchMakerData);
#endif

void create_MPI_type_for_Particles(MPI_Datatype *mpi_particle_type);
void write_grid_to_file(float_kind *grid, double a, char *prefix);

#endif
