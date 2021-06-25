#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include "io_ramses.h"
#include "io_gadget.h"
#include "io_ascii.h"
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

using namespace std;

//======================================
// Simple code to estimate P(k) 
// Hans A. Winther 2015
//======================================

#if defined(_NGP)
const std::string gridassignment = "NGP";
#define powwindow(x) (x)
#elif defined(_CIC)
const std::string gridassignment = "CIC";
#define powwindow(x) ((x)*(x))
#elif defined _TSC
const std::string gridassignment = "TSC";
#define powwindow(x) ((x)*(x)*(x))
#endif


//======================================
// Global container 
//======================================

struct Simulation {
  unsigned long long npart_tot;        // Total number of particles [updated when reading]
  int nfiles;           // Number of input files
  int ngrid;            // FFT grid size
  int nbuffer;          // Size of buffer
  double boxsize;       // Size of the simulation box 
  double *readbuffer;   // Read buffer
  double *pofk;         // P(k) bins
  double *nmodes;       // Number of modes in each pofk bin
  string filebase;      // "/path/to/output_0000X/part_0000X.out" or "/path/to/gadget."
  string pofkoutfile;   // Name of pofk output file
  string datatype;      // RAMSES or GADGET
  fftw_complex *grid;   // FFT grid
  fftw_plan plan;       // FFT execution plan
  bool endianchange;    // finds if the simulation data differs in endianness from the system and changes if need  
} global;

//=================================================
// Nearest-Grid-Point density assignment scheme
//=================================================
void add_to_grid_NGP(double x, double y, double z, int ngrid){
  long int ix, iy, iz;

  ix = (long int)(x * ngrid);
  iy = (long int)(y * ngrid);
  iz = (long int)(z * ngrid);
  
  if(ix >= ngrid) ix -= ngrid;
  if(iy >= ngrid) iy -= ngrid;
  if(iz >= ngrid) iz -= ngrid;

  global.grid[ix + ngrid*(iy + ngrid*iz)][0] += 1.0;
}

//=================================================
// Cloud-In-Cell density assignment scheme
//=================================================
void add_to_grid_CIC(double x, double y, double z, int ngrid){
  long int ix, iy, iz;

  ix = (long int)(x * ngrid);
  iy = (long int)(y * ngrid);
  iz = (long int)(z * ngrid);

  double DX = (x * ngrid) - (double)ix;
  double DY = (y * ngrid) - (double)iy;
  double DZ = (z * ngrid) - (double)iz;

  double TX = 1.0 - DX;
  double TY = 1.0 - DY;
  double TZ = 1.0 - DZ;

  if(ix >= ngrid) ix -= ngrid;
  if(iy >= ngrid) iy -= ngrid;
  if(iz >= ngrid) iz -= ngrid;

  long int ixneigh = ix+1;
  long int iyneigh = iy+1;
  long int izneigh = iz+1;

  if(ixneigh >= ngrid) ixneigh -= ngrid;
  if(iyneigh >= ngrid) iyneigh -= ngrid;
  if(izneigh >= ngrid) izneigh -= ngrid;

  global.grid[ix      + ngrid*(iy      + ngrid*iz     )][0] += TX*TY*TZ;
  global.grid[ix      + ngrid*(iy      + ngrid*izneigh)][0] += TX*TY*DZ;
  global.grid[ix      + ngrid*(iyneigh + ngrid*iz     )][0] += TX*DY*TZ;
  global.grid[ix      + ngrid*(iyneigh + ngrid*izneigh)][0] += TX*DY*DZ;
  
  global.grid[ixneigh + ngrid*(iy      + ngrid*iz     )][0] += DX*TY*TZ;
  global.grid[ixneigh + ngrid*(iy      + ngrid*izneigh)][0] += DX*TY*DZ;
  global.grid[ixneigh + ngrid*(iyneigh + ngrid*iz     )][0] += DX*DY*TZ;
  global.grid[ixneigh + ngrid*(iyneigh + ngrid*izneigh)][0] += DX*DY*DZ;
}



//===================================================
// Triangular-Shaped-Cloud density assignment scheme 
//===================================================

void add_to_grid_TSC(double x, double y, double z, int ngrid){
  long int ix, iy, iz;

  ix = (long int)(x * ngrid);
  iy = (long int)(y * ngrid);
  iz = (long int)(z * ngrid);

  // these are always positive:

  double DX = (x * ngrid) - (double)ix;
  double DY = (y * ngrid) - (double)iy;
  double DZ = (z * ngrid) - (double)iz;
  
  // In each dimension: 
  // N=Next, P=Previous 
  // The nearest grid point gets the weight 0.75 - d^2 
  // Since casting to int rounds down, inbetween two consecutive points, it chooses the former one.
  // The point that follows it gets the weight  0.5*(1.5 - (1-d))^2 = 0.5*(0.5 + d)^2
  // The point that precedes it gets the weight 0.5*(1.5 - (1+d))^2 = 0.5*(0.5 - d)^2 
  // Notice that all is consistent weight-wise  - 
  // i.e. even though the distance gets rounded down, the point that follows the NGP gets more contribution -
  // than the point preceding the NGP - because it is closer.
  // T = THIS, N = NEXT, P = PREVIOUS
  // This scheme sums up to 1.


  double TX = 0.75 - DX*DX;
  double TY = 0.75 - DY*DY;
  double TZ = 0.75 - DZ*DZ;

  double NX = 0.5*(0.5 + DX)*(0.5+DX);
  double NY = 0.5*(0.5 + DY)*(0.5+DY);
  double NZ = 0.5*(0.5 + DZ)*(0.5+DZ);

  double PX = 0.5*(0.5 - DX)*(0.5-DX);
  double PY = 0.5*(0.5 - DY)*(0.5-DY);
  double PZ = 0.5*(0.5 - DZ)*(0.5-DZ);

  // we need to use DX, DY, DZ, NX, NY, NZ, PX, PY, PZ multiplicatively in the final weight assignments = 3**3 possibilities 

  // periodicity for NGPs 
  if(ix >= ngrid) ix -= ngrid;
  if(iy >= ngrid) iy -= ngrid;
  if(iz >= ngrid) iz -= ngrid;

  long int ixneighN = ix+1;  // next 
  long int iyneighN = iy+1;
  long int izneighN = iz+1;

  long int ixneighP = ix-1;  // previous 
  long int iyneighP = iy-1;
  long int izneighP = iz-1;

  // periodicity for points after NGPs

  if(ixneighN >= ngrid) ixneighN -= ngrid;
  if(iyneighN >= ngrid) iyneighN -= ngrid;
  if(izneighN >= ngrid) izneighN -= ngrid;

  // periodicity for points before NGPs

  if(ixneighP < 0) ixneighP += ngrid;
  if(iyneighP < 0) iyneighP += ngrid;
  if(izneighP < 0) izneighP += ngrid;

  // 27 combinations... coming from lowest indices to highest, like this, in the ternary numeral system: 000 001 002 010 011 012 022 100 etc...
  // where 0=P=PREVIOUS, 1=T=MIDDLE(THIS), 2=N=NEXT 
  
  global.grid[ixneighP      + ngrid*(iyneighP      + ngrid*izneighP     )][0] += PX*PY*PZ;   //000
  global.grid[ixneighP      + ngrid*(iyneighP      + ngrid*iz           )][0] += PX*PY*TZ;   //001
  global.grid[ixneighP      + ngrid*(iyneighP      + ngrid*izneighN     )][0] += PX*PY*NZ;   //002

  global.grid[ixneighP      + ngrid*(iy            + ngrid*izneighN     )][0] += PX*TY*PZ;   //010
  global.grid[ixneighP      + ngrid*(iy            + ngrid*iz           )][0] += PX*TY*TZ;   //011
  global.grid[ixneighP      + ngrid*(iy            + ngrid*izneighN     )][0] += PX*TY*NZ;   //012

  global.grid[ixneighP      + ngrid*(iyneighN      + ngrid*izneighP     )][0] += PX*NY*PZ;   //020
  global.grid[ixneighP      + ngrid*(iyneighN      + ngrid*iz           )][0] += PX*NY*TZ;   //021
  global.grid[ixneighP      + ngrid*(iyneighN      + ngrid*izneighN     )][0] += PX*NY*NZ;   //022



  global.grid[ix            + ngrid*(iyneighP      + ngrid*izneighP     )][0] += TX*PY*PZ;   //100
  global.grid[ix            + ngrid*(iyneighP      + ngrid*iz           )][0] += TX*PY*TZ;   //101
  global.grid[ix     	    + ngrid*(iyneighP      + ngrid*izneighN     )][0] += TX*PY*NZ;   //102

  global.grid[ix            + ngrid*(iy            + ngrid*izneighN     )][0] += TX*TY*PZ;   //110
  global.grid[ix            + ngrid*(iy            + ngrid*iz           )][0] += TX*TY*TZ;   //111
  global.grid[ix            + ngrid*(iy            + ngrid*izneighN     )][0] += TX*TY*NZ;   //112

  global.grid[ix            + ngrid*(iyneighN      + ngrid*izneighP     )][0] += TX*NY*PZ;   //120
  global.grid[ix            + ngrid*(iyneighN      + ngrid*iz           )][0] += TX*NY*TZ;   //121
  global.grid[ix            + ngrid*(iyneighN      + ngrid*izneighN     )][0] += TX*NY*NZ;   //122


  global.grid[ixneighN      + ngrid*(iyneighP      + ngrid*izneighP     )][0] += NX*PY*PZ;   //200
  global.grid[ixneighN      + ngrid*(iyneighP      + ngrid*iz           )][0] += NX*PY*TZ;   //201
  global.grid[ixneighN      + ngrid*(iyneighP      + ngrid*izneighN     )][0] += NX*PY*NZ;   //202

  global.grid[ixneighN      + ngrid*(iy            + ngrid*izneighN     )][0] += NX*TY*PZ;   //210
  global.grid[ixneighN      + ngrid*(iy            + ngrid*iz           )][0] += NX*TY*TZ;   //211
  global.grid[ixneighN      + ngrid*(iy            + ngrid*izneighN     )][0] += NX*TY*NZ;   //212

  global.grid[ixneighN      + ngrid*(iyneighN      + ngrid*izneighP     )][0] += NX*NY*PZ;   //220
  global.grid[ixneighN      + ngrid*(iyneighN      + ngrid*iz           )][0] += NX*NY*TZ;   //221
  global.grid[ixneighN      + ngrid*(iyneighN      + ngrid*izneighN     )][0] += NX*NY*NZ;   //222


}

//=================================================
// Process GADGET data
// Format of pos is [x1 y1 z1 z2 y2 z2 ...]
// with x,y,z in [0, BoxSize]
//
// Process RAMSES data
// Format pos [x1 ... xN y1 ... z1 ... ]
// with x,y,z in [0,1]
//=================================================

void process_particle_data(double *pos, int npart, double boxsize = 0.0){
  double x, y, z;
  int ngrid = global.ngrid;
  float *pos_float = (float *) pos;
  

  global.boxsize=boxsize;                      
  #ifdef _UNITS
  global.boxsize/=units;
  #endif


  for(int i = 0; i < npart; i++){
    if(global.datatype.compare("RAMSES") == 0){
      x = pos[i];
      y = pos[npart+i];
      z = pos[2*npart+i];
    } else if (global.datatype.compare("GADGET") == 0){
        if(endianchange){
            swap_endian(pos_float[3*i]);
            swap_endian(pos_float[3*i+1]);
            swap_endian(pos_float[3*i+2]);
            }

      x = double(pos_float[3*i]/boxsize);
      y = double(pos_float[3*i+1]/boxsize);
      z = double(pos_float[3*i+2]/boxsize);
    } else {
      x = pos[3*i]/boxsize;
      y = pos[3*i+1]/boxsize;
      z = pos[3*i+2]/boxsize;
    }


    //======================================
    // Bin particle to the grid
    //======================================
#if defined(_NGP)
    add_to_grid_NGP( x, y, z, ngrid);
#elif defined(_CIC)
    add_to_grid_CIC( x, y, z, ngrid);
#elif defined(_TSC)
    add_to_grid_TSC( x, y, z, ngrid);
#endif
  }

}

//======================================
// Allocate memory and initialize arrays
//======================================
void allocate_memory(bool allocate){
  int ngrid = global.ngrid;

  if(allocate){
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());

    global.grid   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ngrid*ngrid*ngrid);
    global.plan   = fftw_plan_dft_3d(ngrid, ngrid, ngrid, global.grid, global.grid, FFTW_FORWARD, FFTW_ESTIMATE); 
    global.pofk   = new double[ngrid];
    global.nmodes = new double[ngrid];

    // Init P(k) arrays
    for(int i = 0; i < ngrid; i++) 
      global.pofk[i] = global.nmodes[i] = 0.0;

    // Init FFT grid
    int ngridtot = pow3(ngrid);
    for(int i = 0; i < ngridtot; i++)
      global.grid[i][0] = global.grid[i][1] = 0.0;

  } else {
    fftw_free(global.grid);
    delete[] global.pofk;
    delete[] global.nmodes;
  }
}

//======================================
// Window function for the density
// assignment
//======================================
double window(int kx, int ky, int kz){
  double fac, kkx, kky, kkz, w;

  fac = M_PI/double(global.ngrid);
  kkx = kx*fac;
  kky = ky*fac;
  kkz = kz*fac;

  w = 1.0;
  if(kx != 0)
    w *= sin(kkx)/kkx;
  if(ky != 0)
    w *= sin(kky)/kky;
  if(kz != 0)
    w *= sin(kkz)/kkz;

  return powwindow(w);
}

//======================================
// Compute P(k) by FFTing the density field
// Assuming we have binned particles to grid
//======================================
void calculate_pofk(){
  long long ngrid = global.ngrid;
  fftw_complex *grid = global.grid;
  fftw_plan *plan = &global.plan;
  double *pofk    = global.pofk;
  double *nmodes  = global.nmodes;
  double *pofk_t, *nmodes_t;
  int nthreads, id;

  //======================================
  // Make one pofk array for each thread
  //======================================
#pragma omp parallel private(id)
  {
    id = omp_get_thread_num(); 
    if(id == 0) nthreads = omp_get_num_threads();
  }
  pofk_t   = new double[ngrid*nthreads];
  nmodes_t = new double[ngrid*nthreads];
  for(int i = 0; i < ngrid*nthreads; i++)
    pofk_t[i] = nmodes_t[i] = 0.0;

  //======================================
  // Perform FFT
  //=====================================-
  cout << "==> Performing FFT\n" << endl;
  fftw_execute(*plan);

  //======================================
  // Bin up pofk
  //======================================
  cout << "===> Binning up P(k)\n" << endl;
  int nover2 = ngrid/2;
  double fftw_norm_fac = pow( 1.0/double(ngrid*ngrid*ngrid), 2);
  for(int i = 0; i < ngrid; i++){
    int ii = (i < nover2 ? i : i-ngrid);
    for(int j = 0; j < ngrid; j++){
      int jj = (j < nover2 ? j : j-ngrid);
#pragma omp parallel for 
      for(int k = 0; k < ngrid; k++){
        int dind = ngrid * omp_get_thread_num();
        int kk = (k < nover2 ? k : k-ngrid);
        long long ind = i + ngrid*(j + k*ngrid);
        long long int kind = int(sqrt(ii*ii + jj*jj + kk*kk) + 0.5);

        if(kind < ngrid && kind > 0){
          pofk_t[kind+dind] += (grid[ind][0]*grid[ind][0]+grid[ind][1]*grid[ind][1])*fftw_norm_fac / pow2( window(ii,jj,kk) );
          nmodes_t[kind+dind] += 1.0;
        }
      }
    }
  }

  //======================================
  // Sum up pofk from all the threads
  //======================================
  for(int i = 0; i < ngrid; i++){
    for(int j = 0; j < nthreads; j++){
      pofk[i] += pofk_t[i + ngrid*j];
      nmodes[i] += nmodes_t[i + ngrid*j];
    }
  }
  for(int i = 0; i < ngrid; i++){
    if(nmodes[i] > 0){
      pofk[i] /= nmodes[i];
    }
  }

  //======================================
  // Subtract shotnoise
  //======================================
#ifdef _SUBTRACTSHOTNOISE
  for(int i = 0; i < ngrid; i++){
    pofk[i] -= 1.0/double(global.npart_tot);
  }
#endif

  //======================================
  // Clean up
  //======================================
  delete[] pofk_t;
  delete[] nmodes_t;
}
  
//======================================
// Output pofk up to the Nyquist frequency 
// (inaccurate otherwise)
//======================================

void output_pofk(){
  int ngrid      = global.ngrid;
  double *pofk   = global.pofk;
  double *nmodes = global.nmodes;
  double boxsize = global.boxsize;

  ofstream pofkfile(global.pofkoutfile.c_str());
  for(int i = 1; i <= ngrid/2; i++){
    if(nmodes[i]>0){
      pofkfile << (2*i-1)*M_PI/boxsize << " " << double(pofk[i]*boxsize*boxsize*boxsize) << endl;
    }
  } 
}

int main(int argv, char **argc){

  //======================================
  // Initialize parameters
  //======================================
  if(argv < 5){
    cout << "Run as ./pofk /path/output_0000X/part_0000X.out outfilename ngrid nFiles RAMSES/GADGET" << endl;
    exit(1);
  } else {
    global.filebase    = argc[1];
    global.pofkoutfile = argc[2];
    global.ngrid       = atoi(argc[3]);
    global.nfiles      = atoi(argc[4]);
    global.datatype    = argc[5];
    global.npart_tot   = 0;
//  global.nbuffer     = 100000000; // 750 MB
    global.nbuffer     =100000000*10;

    cout << endl;
    cout << "=====================================" << endl;
    cout << "Parameters:                          " << endl;
    cout << "=====================================" << endl;
    cout << "Filebase   = " << global.filebase      << endl; 
    cout << "Outfile    = " << global.pofkoutfile   << endl; 
    cout << "Ngrid      = " << global.ngrid         << endl; 
    cout << "Nfiles     = " << global.nfiles        << endl; 
    cout << "Datatype   = " << global.datatype      << endl; 
    cout << "Assignment = " << gridassignment       << endl;
#pragma omp parallel
    {
      if(omp_get_thread_num() == 0)
        cout << "Nthreads = " << omp_get_num_threads()  << endl;
    }
    cout << "=====================================" << endl;
    cout << endl;
  }

  //======================================
  // Allocate memory
  //======================================
  allocate_memory(true);

  //======================================
  // Read particles from data files
  //======================================
  global.readbuffer = new double[global.nbuffer];
  for(int i = 1; i <= global.nfiles; i++){

    if(global.datatype.compare("RAMSES") == 0){
      cout << "Read " << i << endl;
      read_and_bin_particles_ramses(global.filebase, i  , &global.npart_tot, global.readbuffer, &global.nbuffer);
    } else if(global.datatype.compare("GADGET") == 0) {
      read_and_bin_particles_gadget(global.filebase, i-1, &global.npart_tot, global.readbuffer, &global.nbuffer, global.nfiles);
    } else if(global.datatype.compare("ASCII") == 0) {
      global.nfiles = 1;
      read_and_bin_particles_ascii(global.filebase, 1, &global.npart_tot, global.readbuffer, &global.nbuffer);
    } else {
      cout << "Error: Unknown datatype [" << global.datatype << "]" << endl;
      exit(1);
    }

  }
  delete[] global.readbuffer;

  //======================================
  // Change from rho to delta = rho/rho_b - 1
  //======================================
  cout << "\n=> Normalizing density field\n" << endl;
  
  long long ngrid = global.ngrid;  
  double fac = double(ngrid*ngrid*ngrid)/double(global.npart_tot);  
  double avg = 0.0, maxdens = 0.0;

  for(unsigned long long i=0; i<(unsigned long long)(ngrid*ngrid*ngrid);i++){
  	global.grid[i][0]/=fac;
  	avg += global.grid[i][0];	    
	}
  
  avg /= double(ngrid*ngrid*ngrid);

  for(unsigned long long i=0; i<(unsigned long long)(ngrid*ngrid*ngrid);i++){
        if(global.grid[i][0] > maxdens) maxdens = global.grid[i][0];
	global.grid[i][0] = global.grid[i][0]/avg - 1.0;
	}

  printf("Average density rho = %f  Maxdens = %f\n", avg, maxdens);

  //======================================
  // Compute and output P(k)
  //======================================
  calculate_pofk();
  output_pofk();
  
  //======================================
  // Clean up
  //======================================
  allocate_memory(false);
}

