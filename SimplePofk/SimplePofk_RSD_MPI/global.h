#define FILETYPE_OWNFORMAT 0
#define FILETYPE_GADGET    1
#define FILETYPE_ASCII     2

#define LINEAR_SPACING     0
#define LOG_SPACING        1

#define FALSE              0
#define TRUE               1

//====================================
// This is position_unit / (Mpc/h)
// If positions are in Mpc/h then 1.0
// If kpc/h then use 1000.0
//====================================
#define MPC_UNIT 1.0

//==============================================
// Types used for the arrays. If one wants to use 
// float to save memory then fftw -> fftwf 
// in the routines below must be changed
//==============================================
typedef double float_kind;
typedef fftw_complex complex_kind;

// 0 = x, 1 = y, 2 = z
int line_of_sight_direction = 0;

double OmegaM;
double OmegaLambda;
double redshift;

//==============================================
// MPI variables
//==============================================
struct MPIVariables{
  ptrdiff_t Local_nx;     
  ptrdiff_t Total_size;   
  ptrdiff_t alloc_local;  
  ptrdiff_t alloc_slice;  
  ptrdiff_t Local_x_start;
  ptrdiff_t last_slice;
  int *Local_nx_table;
  int NTask;
  int ThisTask;
  int LeftTask, RightTask;
  int *Slab_to_task;
  int ierr;
  MPI_Status status; 
} mpi;

//==============================================
// Particle struct
//==============================================
struct Particle {
  float_kind Pos[3];
};

struct Particles {
  ptrdiff_t NumPart;
  ptrdiff_t NumPartTot;
  struct Particle *P;
  double buffer;
};

//==============================================
// Global parameters
//==============================================
struct Global {
  double buffer;
  double Box;
  double OmegaM;
  double OmegaLambda;
  double redshift;
  int los;
} global;

//==============================================
// Settings for P(k) estimation
//==============================================
struct PowerSpectrum{
  int nbins;               // Number of bins
  int bintype;             // The bintype (linear / logspaced)
  int subtract_shotnoise;  // Subtract shotnoise or not
  double kmin;             // Integer wavenumber (0 <= kmin <= Nmesh)
  double kmax;             // Integer wavenumber (kmin < kmax <= Nmesh)
  char outname[100];       // Output filename
} pofksettings;

//==============================================
// All methods below
//==============================================
void   initialize_ffts(int Nmesh);
void   initialize_particles(char *filebase, struct Particles *p, int Nmesh);
void   initialize_particles_gadget(char *filebase,  struct Particles *P, int Nmesh);
void   initialize_particles_ascii(char *filebase, struct Particles *p, int Nmesh);
void   ParticlesToGrid(struct Particles *p, float_kind *density, int Nmesh);
int    pofk_bin_index(double kmag, double kmin, double kmax, int nbins, int bintype);
void   compute_2d_pofk(struct Particles *p, complex_kind *dens_k, int Nmesh);
double k_from_index(int index, double kmin, double kmax, int nbins, int bintype);

#ifdef _CENTRALOBSERVER

void make_source_centralobserver(float_kind *density, int Nmesh);
void add_source(complex_kind *density_source, complex_kind *density_ij, int Nmesh, int ind_1, int ind_2);
void compute_density_ij(float_kind *density, float_kind *density_xixj, int Nmesh, int ind_1, int ind_2);

#endif
