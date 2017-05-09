#include <fftw3.h>
#include <fftw3-mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "io_gadget.h"

//==================================================================\\
//                                                                  \\
// Compute the RSD matter power-spectrum multipoles                 \\
// Written by Hans Winther                                          \\
//                                                                  \\
// Each CPU has all particles for which                             \\
//                                                                  \\
// mpi.Local_x_start/Nmesh < X/Box <= mpi.Local_nx/Nmesh            \\
//                                                                  \\
// where Box is the boxsize. If this is not in Mpc/h one            \\
// needs to set MPC_UNIT appropriately                              \\
//                                                                  \\
// The function initialize_particles() needs to be implemented      \\
// for the particular file-format one uses. GADGET1 format          \\
// implemented so far.                                              \\
//                                                                  \\
// The output is k (h/Mpc), P(k) (Mpc/h)^3                          \\
// NB: This assumes the BoxSize in the input-files is in [Mpc/h]    \\
//                                                                  \\
// The other variables one needs to set are                         \\
// * Box (the boxsize, i.e. 0 <= X,Y,Z < Box)                       \\
// * Nmesh (the 1D size of the grid to compute P(k) on)             \\
// * NumPartTot (the total number of particles)                     \\
// * buffer (allocate x times the mean number of particles per cpu) \\
//                                                                  \\
// Modify the pofksettings variables to change the binning          \\
// if wanted.                                                       \\
//                                                                  \\
// The code should hopefully work for large particle-numbers,       \\
// as long as we run with enough processes, but it has not been     \\
// tested for this case.                                            \\
//                                                                  \\
//==================================================================\\

//==============================================
// Set up the sizes for arrays and fft routines
//==============================================
void initialize_ffts(int Nmesh) {
  int *Slab_to_task_local;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.NTask);
  fftw_mpi_init();

  mpi.alloc_local = fftw_mpi_local_size_3d(Nmesh, Nmesh, Nmesh/2+1, MPI_COMM_WORLD, &mpi.Local_nx, &mpi.Local_x_start);
  mpi.Local_nx_table = malloc(sizeof(int) * mpi.NTask);
  MPI_Allgather(&mpi.Local_nx, 1, MPI_INT, mpi.Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);
  if(mpi.ThisTask == 0) {
    printf("\nLocal nodes per Task:\n---------------------\n");
    for(int i = 0; i < mpi.NTask; i++) printf("Task = %d: mpi.Local_nx = %d\n", i, mpi.Local_nx_table[i]);
    printf("---------------------\n\n");
    fflush(stdout);
  }

  //==============================================
  // Set the neighbouring tasks
  //==============================================
  if (mpi.Local_nx == 0) {
    mpi.LeftTask = MPI_PROC_NULL;
    mpi.RightTask = MPI_PROC_NULL;
  } else {
    mpi.LeftTask = mpi.ThisTask;
    do {
      mpi.LeftTask--;
      if(mpi.LeftTask < 0) mpi.LeftTask = mpi.NTask - 1;
    } while(mpi.Local_nx_table[mpi.LeftTask] == 0);

    mpi.RightTask = mpi.ThisTask;
    do {
      mpi.RightTask++;
      if(mpi.RightTask >= mpi.NTask) mpi.RightTask = 0;
    } while(mpi.Local_nx_table[mpi.RightTask] == 0);
  }

  //==============================================
  // Let each processor know which parts of the fourier grids they all have
  //==============================================
  mpi.Slab_to_task       = malloc(sizeof(int) * Nmesh);
  Slab_to_task_local = malloc(sizeof(int) * Nmesh);

  for(int i = 0; i < Nmesh; i++) Slab_to_task_local[i] = 0;
  for(int i = 0; i < mpi.Local_nx; i++) Slab_to_task_local[mpi.Local_x_start + i] = mpi.ThisTask;

  MPI_Allreduce(Slab_to_task_local, mpi.Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  //==============================================
  // Add an additional plane
  //==============================================
  mpi.alloc_slice = Nmesh*(Nmesh/2+1);
  mpi.last_slice = mpi.Local_nx * mpi.alloc_slice;
  mpi.Total_size = mpi.alloc_local + mpi.alloc_slice;

  free(Slab_to_task_local);
  return;
}

//===========================================================
// Read particles from file. Assume positions are in [0,Box]
// After this routine is done Part.P should contain all the 
// particles with x-positions satisfying:
//
// mpi.Local_x_start/Nmesh < X/Box <= mpi.Local_nx/Nmesh 
//
// and NumPart should contain the number of particles on
// the current CPU and NumPartTot the total over all CPUs.
// 
// See initialize_particles_ascii or _gadget for how to 
// write this routine
//===========================================================

void initialize_particles(char *filebase, struct Particles *p, int Nmesh){
  double scaleBox = (double)Nmesh / p->Box;
  struct Particle *part = p->P;

  // This method needs to be written by the user...
  printf("Method initialize_particles is not implented...\n");
  MPI_Abort(MPI_COMM_WORLD,1);
  exit(1);

  // Set the total number of particles across all CPUs (used to normalize the density-field)
  p->NumPartTot = 10000;

  // Estimate for the number of particles
  p->NumPart = (unsigned int)( (double) p->NumPartTot / (double) mpi.NTask);
  
  // Allocate memory for particles
  ptrdiff_t NPartAllocated = (ptrdiff_t)(p->NumPart * p->buffer);
  part = malloc((ptrdiff_t)( NPartAllocated * sizeof(struct Particle)));
  

  // Open file(s)
  // ...
  // ...
  // ...

  // Read files and add particles that are within the slice
  ptrdiff_t NumPart_read = 0;
  for(;;){
    double X = 0.0, Y = 0.0, Z = 0.0;
    
    // Read X, Y, Z coordinates of current particle
    // ...
    // ...
    // ...

    // We only process particles that are in the slice belonging to this CPU
    int IX = (int)(X * scaleBox) - (int)(mpi.Local_x_start);
    if( IX >= mpi.Local_nx || IX < 0) continue;

    // Assign particles
    part[NumPart_read].Pos[0] = X;
    part[NumPart_read].Pos[1] = Y;
    part[NumPart_read].Pos[2] = Z;
    NumPart_read++;

    // Sanity check
    if(NumPart_read >= NPartAllocated){
      printf("Error in reading particles. Increase buffer\n");
      MPI_Abort(MPI_COMM_WORLD,1);
      exit(1);
    }
  }

  // The acctual number of particles we have in the slice
  p->NumPart = NumPart_read;
}

//===========================================================
// Read gadget files
//===========================================================

// Book-keeping variable to keep track of how many particles we have read that
// belongs to this task. Updated in process_particle_data
unsigned int nread_gadget = 0;

void initialize_particles_gadget(char *filebase, struct Particles *p, int Nmesh){
  int nfiles = 0;

  // Read header of files and get total number of particles
  int ifile = 0;
  unsigned int npartmax = 0;
  p->NumPartTot = 0;
  for(;;){
    char filename[100];
    sprintf(filename, "%s%i", filebase, ifile);
    FILE *fp = fopen(filename, "r");
    if(fp == NULL){
      printf("Error: cannot find file [%s]\n", filename);
      MPI_Abort(MPI_COMM_WORLD,1);
      exit(1);
    }
    read_gadget_header(fp);
    if(ifile == 0) nfiles = header.num_files;
    p->NumPartTot += (ptrdiff_t)(header.npart[1]);
    if(header.npart[1] > npartmax) npartmax = header.npart[1];
    fclose(fp);
    if(ifile++ == nfiles - 1) break;
  }
  
  // Make a read-buffer. Must be larger than 3 * npart_max_per_file
  ptrdiff_t nbuffer = 3 * npartmax;
#ifndef _READVELOCITY
  // For real-space P(k) we don't need velocities
  float *readbuffer = malloc(sizeof(float) * nbuffer);
#else
  float *readbuffer = malloc(sizeof(float) * nbuffer * 2);
#endif

  // Allocate memory for particles
  p->NumPart = (ptrdiff_t)( (double) p->NumPartTot / (double) mpi.NTask);
  p->P = malloc((ptrdiff_t)(p->NumPart * p->buffer) * sizeof(struct Particle));

  nread_gadget = 0;
  for(int i = 0; i < nfiles; i++){
    read_and_bin_particles_gadget(filebase, i, &p->NumPartTot, readbuffer, &nbuffer, p, Nmesh);
  }
  p->NumPart = nread_gadget;

  // Show how many particles we have per task
  if(mpi.ThisTask == 0) printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Done reading. ThisTask: [%i] NumPart: [%td]   NumPartTot: [%td]\n", mpi.ThisTask, p->NumPart, p->NumPartTot);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if(mpi.ThisTask == 0) printf("\n");
  
  free(readbuffer);
}

//===========================================================
// Process a single gadget file
//===========================================================
void process_particle_data(double *pos, double *vel, ptrdiff_t npart, double boxsize, struct Particles *p, int Nmesh){
  double scaleBox = (double) Nmesh / boxsize;
  float_kind X, Y, Z;
  struct Particle *part = p->P; 
  float *pos_float = (float *) pos;
#ifdef _READVELOCITY
  float *vel_float = (float *) vel;
#endif

  // Set boxsize from the value found in the files
  p->Box = boxsize;
 
  // Process particles
  for(unsigned int i = 0; i < npart; i++){
    X = (float_kind)(pos_float[3*i]   );
    Y = (float_kind)(pos_float[3*i+1] );
    Z = (float_kind)(pos_float[3*i+2] );
    
    // We only process particles that are in the slice belonging to this CPU
    int IX = (int)(X * scaleBox) - (int)(mpi.Local_x_start);
    if( IX >= mpi.Local_nx || IX < 0) continue;
  
    // Store particle
    part[nread_gadget].Pos[0] = X;
    part[nread_gadget].Pos[1] = Y;
    part[nread_gadget].Pos[2] = Z;
    nread_gadget++;
  }
}

//===========================================================
// Example: read a single ascii-file with format
// npart_tot
// x y z vx vy vz
// x y z vx vy vz
// ...
//===========================================================
void initialize_particles_ascii(char *filebase, struct Particles *p, int Nmesh){
  double scaleBox = (double)Nmesh / p->Box;
  struct Particle *part = p->P;

  // Open file
  FILE* fp = fopen(filebase, "r");
  if(fp == NULL){
    printf("Error: cannot open file [%s]\n", filebase);
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  // Read total number of particles
  fscanf(fp,"%td\n", &p->NumPartTot);
  if(mpi.ThisTask == 0) printf("The particle file has [%td] particles\n", p->NumPartTot);

  // Estimate for the number of particles
  p->NumPart = (unsigned int)(p->NumPartTot / (double) mpi.NTask);
  
  // Allocate memory for particles
  ptrdiff_t NPartAllocated = (ptrdiff_t)(p->NumPart * p->buffer);
  part = malloc((ptrdiff_t)( NPartAllocated * sizeof(struct Particle)));

  // Each CPU reads all particles and keeps only the ones it need
  ptrdiff_t NumPart_read = 0;
  for(ptrdiff_t i = 0; i < p->NumPartTot; i++){
    double v[3], X, Y, Z;
    
    // Format of line is [X Y Z Vx Vy Vz]
    fscanf(fp,"%lf %lf %lf %lf %lf %lf\n", &X, &Y, &Z, &v[0], &v[1], &v[2]); 

    // We only process particles that are in the slice belonging to this CPU
    int IX = (int)(X * scaleBox) - (int)(mpi.Local_x_start);
    if( IX >= mpi.Local_nx || IX < 0) continue;
    
    // Assign particles
    part[NumPart_read].Pos[0] = X;
    part[NumPart_read].Pos[1] = Y;
    part[NumPart_read].Pos[2] = Z;
    NumPart_read++;

    // Sanity check
    if(NumPart_read >= NPartAllocated){
      printf("Error in reading particles. Increase buffer\n");
      MPI_Abort(MPI_COMM_WORLD,1);
      exit(1);
    }
  }
  fclose(fp);
  p->NumPart = NumPart_read;

  printf("Cpu [%i] has [%td] particles\n", mpi.ThisTask, p->NumPart);
}

//==============================
// Does Cloud-in-Cell assignment
//==============================
void ParticlesToGrid(struct Particles *p, float_kind *density, int Nmesh) {
  ptrdiff_t npart = p->NumPart;
  ptrdiff_t nparttot = p->NumPartTot;
  unsigned int IX, IY, IZ;
  unsigned int IXneigh, IYneigh, IZneigh;
  double X, Y, Z;
  double TX, TY, TZ;
  double DX, DY, DZ;
  double scaleBox = (double)Nmesh / p->Box;
  double WPAR = pow((double)Nmesh,3) / (double) nparttot;

  // Initialize density to -1
  for(ptrdiff_t i = 0; i < 2 * mpi.Total_size; i++) 
    density[i] = -1.0;
  
  for(ptrdiff_t i = 0; i < npart; i++) {

    // Scale positions to be in [0, Nmesh]
    X = p->P[i].Pos[0] * scaleBox;
    Y = p->P[i].Pos[1] * scaleBox;
    Z = p->P[i].Pos[2] * scaleBox;

    // Grid-index for cell containing particle
    IX = (unsigned int)X;
    IY = (unsigned int)Y;
    IZ = (unsigned int)Z;

    // Coordinate distances to center of cell
    DX = X-(double)IX;
    DY = Y-(double)IY;
    DZ = Z-(double)IZ;
   
    // CIC weights
    TX = 1.0 - DX;
    TY = 1.0 - DY;
    TZ = 1.0 - DZ;
    DY *= WPAR;
    TY *= WPAR;

    // Periodic BC
    IX -= mpi.Local_x_start;
    if(IY >= (unsigned int)Nmesh) IY = 0;
    if(IZ >= (unsigned int)Nmesh) IZ = 0;

    // Neighbor gridindex
    // No check for x as we have an additional slice on the right
    IXneigh = IX + 1;
    IYneigh = IY + 1;
    IZneigh = IZ + 1;
    if(IYneigh >= (unsigned int)Nmesh) IYneigh = 0;
    if(IZneigh >= (unsigned int)Nmesh) IZneigh = 0;

    //====================================================================================
    // Assign density to the 8 cells containing the particle cloud
    //====================================================================================
    ptrdiff_t index = (ptrdiff_t)(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ;
    density[index] += TX*TY*TZ;
    
    index = (ptrdiff_t)(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh;
    density[index] += TX*TY*DZ;

    index = (ptrdiff_t)(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ;
    density[index] += TX*DY*TZ;
    
    index = (ptrdiff_t)(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh;
    density[index] += TX*DY*DZ;
    
    index = (ptrdiff_t)(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ;
    density[index] += DX*TY*TZ;
    
    index = (ptrdiff_t)(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh;
    density[index] += DX*TY*DZ;
    
    index = (ptrdiff_t)(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ;
    density[index] += DX*DY*TZ;
    
    index = (ptrdiff_t)(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh;
    density[index] += DX*DY*DZ;
  }

  //====================================================================================
  // Copy across the extra slice from the task on the left and add it to the leftmost slice
  // of the task on the right. Skip over tasks without any slices.
  //====================================================================================
  float_kind * temp_density = (float_kind *)calloc(2*mpi.alloc_slice,sizeof(float_kind));
  mpi.ierr = MPI_Sendrecv(&(density[2*mpi.last_slice]),2*mpi.alloc_slice*sizeof(float_kind),MPI_BYTE,mpi.RightTask,0,
      &(temp_density[0]),2*mpi.alloc_slice*sizeof(float_kind),MPI_BYTE,mpi.LeftTask,0,MPI_COMM_WORLD,&mpi.status);
  if (p->NumPart != 0) {
    for (ptrdiff_t i = 0; i < 2 * mpi.alloc_slice; i++) density[i] += (temp_density[i] + 1.0);
  }
  free(temp_density);
}

//===========================================================================
// Defines the binning for the power-spectrum estimation
//===========================================================================
int pofk_bin_index(double kmag, double kmin, double kmax, int nbins, int bintype){

  if(bintype == LINEAR_SPACING) {
    // Linear bins
    int index = (int)( (kmag - kmin) / (kmax - kmin) * nbins + 0.5);
    return index;
  } else if(bintype == LOG_SPACING) {
    // Logarithmic bins
    if(kmag <= 0.0) return -1;
    int index = (int)( log(kmag/kmin)/log(kmax/kmin) * nbins + 0.5);
    return index;
  } else {
    printf("Error: unknown bintype\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  return 0;
}
double k_from_index(int index, double kmin, double kmax, int nbins, int bintype){

  if(bintype == LINEAR_SPACING) {
    double kmag = kmin + (kmax - kmin)/(double)(nbins) * index;
    return kmag;
  } else if(bintype == LOG_SPACING) {
    double kmag = exp(log(kmin) + log(kmax/kmin)/(double)(nbins) * index);
    return kmag;
  } else {
    printf("Error: unknown bintype\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  return 0.0;
}

//=======================================================
// Compute the power-spectrum
//=======================================================
void compute_pofk(struct Particles *p, complex_kind *dens_k, int Nmesh){
  ptrdiff_t coord;
  double pofk, kmag, grid_corr_x, grid_corr_y, grid_corr_z, grid_corr = 1.0;
  double PI = acos(-1.0);
  int d[3], nk;

  //=======================================================
  // Define the binning. Bintype: 0 = Linear; 1 = Logspaced
  // Every integer k-mode: kmin = 0, kmax = Nmesh/2, 
  // nbins = Nmesh/2+1 and bintype = 0
  //
  // Here we fetch the settings from the pofk struct
  //=======================================================
  const int nbins              = pofksettings.nbins;
  const int bintype            = pofksettings.bintype;
  const int subtract_shotnoise = pofksettings.subtract_shotnoise;
  const double kmin            = pofksettings.kmin;
  const double kmax            = pofksettings.kmax;
  char *filename               = pofksettings.outname;

  double *pofk_bin     = malloc(sizeof(double) * nbins);
  double *pofk_bin_all = malloc(sizeof(double) * nbins);
  double *n_bin        = malloc(sizeof(double) * nbins);
  double *n_bin_all    = malloc(sizeof(double) * nbins);
  double *k_bin        = malloc(sizeof(double) * nbins);
  double *k_bin_all    = malloc(sizeof(double) * nbins);

  // FFT normalization factor for |density(k)|^2
  double fac = 1.0/(double) (Nmesh * Nmesh * Nmesh);
  fac = fac*fac;

  for(int i = 0; i < nbins; i++){
    pofk_bin[i] = pofk_bin_all[i] = 0.0;
    n_bin[i]    = n_bin_all[i]    = 0.0;
    k_bin[i]    = k_bin_all[i]    = 0.0;
  }

  // Loop over all modes and add up P(k). We use the symmetry F[i,j,k] = F^*[-i-j,-k] to 
  // get the negative k modes that FFTW does not store
  for (int i = 0; i < mpi.Local_nx; i++) {
    int iglobal = i + mpi.Local_x_start;

    // |kx| component
    d[0] = iglobal > Nmesh/2 ? Nmesh - iglobal : iglobal;

#ifdef _DECONVOLVEWINDOWFUNCTION
    // Window-function (Sqrt[W]) for kx
    grid_corr_x = d[0] == 0 ? 1.0 : sin((PI*d[0])/(double)Nmesh)/((PI*d[0])/(double)Nmesh);
#else
    grid_corr_x = 1.0;
#endif

    for (int j = 0 ; j < (unsigned int)(Nmesh/2+1); j++) {

      // |ky| component
      d[1] = j;

#ifdef _DECONVOLVEWINDOWFUNCTION
      // Window-function (Sqrt[W]) for ky
      grid_corr_y = d[1] == 0 ? 1.0 : sin((PI*d[1])/(double)Nmesh)/((PI*d[1])/(double)Nmesh);
#else
      grid_corr_y = 1.0;
#endif

      //============================
      // Do the k = 0 mode
      //============================
      d[2] = 0;
      grid_corr_z = 1.0;
      kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] ) + 1e-100;
      nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
      if(nk >= 0 && nk < nbins){
        coord = (ptrdiff_t)(i*Nmesh+j)*(Nmesh/2+1) + d[2];
        grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fac;
        pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
        pofk_bin[nk] += 1.0*pofk;
        k_bin[nk]    += 1.0*kmag;
        n_bin[nk]    += 1.0;

        // Mirror around y-axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          coord = (ptrdiff_t)(i*Nmesh+(Nmesh-j))*(Nmesh/2+1) + d[2];
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk_bin[nk] += 1.0*pofk;
          k_bin[nk]    += 1.0*kmag;
          n_bin[nk]    += 1.0;
        }
      }
      
      //============================
      // Do the k = N/2 mode
      //============================
      d[2] = Nmesh/2;
#ifdef _DECONVOLVEWINDOWFUNCTION
      grid_corr_z = 2.0 / PI;
#else
      grid_corr_z = 1.0;
#endif
      kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
      nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
      if(nk >= 0 && nk < nbins){
        coord = (ptrdiff_t)(i*Nmesh+j)*(Nmesh/2+1) + d[2];
        grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fac;
        pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
        pofk_bin[nk] += 1.0*pofk;
        k_bin[nk]    += 1.0*kmag;
        n_bin[nk]    += 1.0;

        // Mirror around y-axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          coord = (ptrdiff_t)(i*Nmesh+(Nmesh-j))*(Nmesh/2+1) + d[2];
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk_bin[nk] += 1.0*pofk;
          k_bin[nk]    += 1.0*kmag;
          n_bin[nk]    += 1.0;
        }
      }

      for (int k = 1; k < (unsigned int) (Nmesh/2); k++) {

        // |kz| component
        d[2] = k;

#ifdef _DECONVOLVEWINDOWFUNCTION
        // Window-function (Sqrt[W]) for kz
        grid_corr_z = sin((PI*d[2])/(double)Nmesh)/((PI*d[2])/(double)Nmesh);
#else
        grid_corr_z = 1.0;
#endif

        //============================
        // Do the general mode
        //============================
        kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
        nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
        if(nk >= 0 && nk < nbins){
          coord = (ptrdiff_t)(i*Nmesh+j)*(Nmesh/2+1) + d[2];
          grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fac;
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk_bin[nk] += 2.0*pofk;
          k_bin[nk]    += 2.0*kmag;
          n_bin[nk]    += 2.0;

          // Mirror around y-axis
          if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
            coord = (ptrdiff_t)(i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + d[2];
            pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
            pofk_bin[nk] += 2.0*pofk;
            k_bin[nk]    += 2.0*kmag;
            n_bin[nk]    += 2.0;
          }
        }
      }
    }
  }

  // Communicate
  MPI_Allreduce(pofk_bin, pofk_bin_all, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(n_bin,    n_bin_all,    nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(k_bin,    k_bin_all,    nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  // Normalize and subtract shotnoise
  for(int i = 0; i < nbins; i++){
    if(n_bin_all[i] > 0){
      pofk_bin[i] = (pofk_bin_all[i] / n_bin_all[i]);
      if(subtract_shotnoise) pofk_bin[i] -= 1.0 / (double) p->NumPartTot;
      k_bin[i]    = (k_bin_all[i]    / n_bin_all[i]);
    }
  }

  // Write output to file
  if(mpi.ThisTask == 0){
    
    // Make filename
    printf("Writing power-spectrum to file: [%s]\n", filename);
    FILE *fp = fopen(filename,"w");
    fprintf(fp,"#     k  (h/Mpc)            P(k)   (Mpc/h)^3   BoxSize: [%f] Mpc/h\n", p->Box);
    for(int i = 1; i < nbins; i++){
      if(n_bin_all[i] > 0){
        double k_of_bin = k_from_index(i, kmin, kmax, nbins, bintype);
        fprintf(fp, "%15.10f   %15.10f\n", k_of_bin * 2.0 * acos(-1.0) / p->Box, pofk_bin[i] * pow(p->Box, 3));
      
        // Alternative: output mean k-value in the bin
        // fprintf(fp, "%14.7f   %14.7f\n", k_bin[i], pofk_bin[i]);
      }
    }
    fclose(fp);
  }

  // Free memory
  free(pofk_bin);
  free(pofk_bin_all);
  free(n_bin);
  free(n_bin_all);
  free(k_bin);
  free(k_bin_all);
}

int main(int argc, char **argv){ 
  struct Particles Part;
    
  if(argc < 4){
    printf("Error: run as ./simplepofk_mpi filebase outfilename ngrid\n");
    printf("Filebase: /path/gadget. for GADGET (only one implemented so far)\n");
    exit(1);
  }

  //====================================================
  // Set the parameters of the run
  //====================================================

  // Set grid-size
  int Nmesh = atoi(argv[3]); 

  // Particles
  char *filebase = argv[1];
  int filetype   = FILETYPE_GADGET;
  Part.buffer    = 2.0; // We allocate x times more memory for the particles than the average

  // Boxsize (used to scale particle positions to [0,1])
  // For GADGET this is taken to be the value in the gadget-header
  Part.Box = 1.0;
  
  // Standard binning of P(k) settings
  char *outputfilename  = argv[2];
  pofksettings.nbins    = Nmesh/2;
  pofksettings.bintype  = LINEAR_SPACING;
  pofksettings.subtract_shotnoise = TRUE;
  pofksettings.kmin     = 0.0;
  pofksettings.kmax     = Nmesh/2;
  sprintf(pofksettings.outname,"%s", outputfilename);

  //====================================================
  
  // Initialize MPI and FFTW
  MPI_Init(&argc, &argv); 

  if(mpi.ThisTask == 0){
    printf("Nmesh:    [%i]\nFilebase: [%s]\nOutput:   [%s]\n", Nmesh, filebase, outputfilename);
  }

  initialize_ffts(Nmesh);

  // Allocate memory and read in particles
  if(mpi.ThisTask == 0) printf("Reading particles from file...\n");
  if(filetype == FILETYPE_OWNFORMAT){
    initialize_particles(filebase, &Part, Nmesh);
  } else if(filetype == FILETYPE_GADGET){
    initialize_particles_gadget(filebase, &Part, Nmesh);
  } else if(filetype == FILETYPE_ASCII){
    initialize_particles_ascii(filebase, &Part, Nmesh);
  } else {
    printf("Filetype not recognized. Exit!\n");
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  // Make density array
  float_kind *density = malloc(2 * mpi.Total_size * sizeof(float_kind));
  complex_kind *density_k = (complex_kind *) density;

  // Bin particles to grid
  if(mpi.ThisTask == 0) printf("Binning particles to grid...\n");
  ParticlesToGrid(&Part, density, Nmesh);

  // Fourier transform the density field
  if(mpi.ThisTask == 0) printf("Fourier transforming density field...\n");
  fftw_plan plan = fftw_mpi_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh, density, density_k, MPI_COMM_WORLD, FFTW_ESTIMATE); 
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // Bin up P(k)
  if(mpi.ThisTask == 0) printf("Binning up P(k)...\n");
  compute_pofk(&Part, density_k, Nmesh);

  // Clean up
  free(density);
  free(Part.P);
  MPI_Finalize();
}

