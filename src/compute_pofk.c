#include "vars.h"
#include "proto.h"
#include "msg.h"
#include "timer.h"
#define LINEAR_SPACING 0
#define LOG_SPACING    1
#define XAXIS          0
#define YAXIS          1
#define ZAXIS          2

//==========================================================================//
//                                                                          //
//  MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017           //
//                                                                          //
//  This files contains methods to compute P(k) and RSD multiplole P_ell(k) //
//                                                                          //
//==========================================================================//

void compute_power_spectrum(complex_kind *dens_k, double a, char *label);
void adjust_pofk_parameters(int *nbins, int *bintype, int *subtract_shotnoise, double *kmin, double *kmax);
void compute_RSD_powerspectrum(double A, int dDdy_set_in_particles);
void PtoMesh_RSD(float_kind *dens, int axis, double A);
void bin_up_RSD_power_spectrum(complex_kind *dens_k, struct RSD_Pofk_Data *pofk_data);
int pofk_bin_index(double kmag, double kmin, double kmax, int nbins, int bintype);
double k_from_index(int index, double kmin, double kmax, int nbins, int bintype);

//===========================================================================
// Defines the binning for the power-spectrum estimation
//===========================================================================
int pofk_bin_index(double kmag, double kmin, double kmax, int nbins, int bintype){

  // Linear bins
  if(bintype == LINEAR_SPACING) {
    int index = (int)( (kmag - kmin) / (kmax - kmin) * nbins + 0.5);
    return index;
  }

  // Logarithmic bins
  if(bintype == LOG_SPACING) {
    if(kmag <= 0.0) return -1;
    int index = (int)( log(kmag/kmin)/log(kmax/kmin) * nbins + 0.5);
    return index;
  }

  return 0;
}

//===========================================================================
// Compute k in h/Mpc from it's index
//===========================================================================
double k_from_index(int index, double kmin, double kmax, int nbins, int bintype){

  if(bintype == LINEAR_SPACING) {
    double kmag = kmin + (kmax - kmin)/(double)(nbins) * index;
    return kmag * 2.0 * M_PI / Box;
  }
  
  if(bintype == LOG_SPACING) {
    double kmag = exp(log(kmin) + log(kmax/kmin)/(double)(nbins) * index);
    return kmag * 2.0 * M_PI / Box;
  }

  return 0.0;
}

//===========================================================================
// Compute P(k,a) = <|density(k,a)|^2>
// Assumes dens_k has the fourier transform of the density field
// The scale-factor is just use to make the output-name
//===========================================================================
void compute_power_spectrum(complex_kind *dens_k, double a, char *label){
  timer_start(_PofkComputation);
  unsigned int coord;
  double pofk, kmag, grid_corr_x, grid_corr_y, grid_corr_z, grid_corr = 1.0;
  int d[3], nk;

  //=======================================================
  // Define the binning. Bintype: 0 = Linear; 1 = Logspaced
  // Every integer k-mode: kmin = 0, kmax = Nmesh/2, 
  // nbins = Nmesh/2+1 and bintype = 0
  //=======================================================
  
#ifdef COMPUTE_POFK

  // Get parameters from parameterfile
  int nbins              = pofk_nbins;                     // Number of bins in k
  int bintype            = pofk_bintype;;                  // Linear [0] Logarithmic [1]
  int subtract_shotnoise = pofk_subtract_shotnoise;        // Shot-noise subtraction
  double kmin            = pofk_kmin * Box / (2.0 * M_PI); // Min 'integer' wave-number
  double kmax            = pofk_kmax * Box / (2.0 * M_PI); // Max 'integer' wave-number

#else

  // Fiducial parameters [Nmesh/2, LINEAR_SPACING, 1, 0.0, Nmesh/2]
  int nbins              = Nmesh/2;
  int bintype            = LINEAR_SPACING;
  int subtract_shotnoise = 1;
  double kmin            = 0.0;
  double kmax            = Nmesh/2.0;

#endif
 
  // Sanity checks. Adjust parameters to sane values
  adjust_pofk_parameters(&nbins, &bintype, &subtract_shotnoise, &kmin, &kmax);

  double *pofk_bin     = my_malloc(sizeof(double) * nbins);
  double *pofk_bin_all = my_malloc(sizeof(double) * nbins);
  double *n_bin        = my_malloc(sizeof(double) * nbins);
  double *n_bin_all    = my_malloc(sizeof(double) * nbins);
  double *k_bin        = my_malloc(sizeof(double) * nbins);
  double *k_bin_all    = my_malloc(sizeof(double) * nbins);

  for(int i = 0; i < nbins; i++){
    pofk_bin[i] = pofk_bin_all[i] = 0.0;
    n_bin[i]    = n_bin_all[i]    = 0.0;
    k_bin[i]    = k_bin_all[i]    = 0.0;
  }

  // FFT normalization factor for |density(k)|^2
  double fftw_norm_fac = 1.0 / pow( (double) Nmesh, 6);

  // Loop over all modes and bin up P(k)
  for (int i = 0; i < Local_nx; i++) {
    int iglobal = i + Local_x_start;

    // |kx| component
    d[0] = iglobal > Nmesh/2 ? Nmesh - iglobal : iglobal;

    // Window-function (Sqrt[W]) for kx
    grid_corr_x = d[0] == 0 ? 1.0 : sin((PI*d[0])/(double)Nmesh)/((PI*d[0])/(double)Nmesh);

    for (int j = 0 ; j < (unsigned int)(Nmesh/2+1); j++) {

      // |ky| component
      d[1] = j;

      // Window-function (Sqrt[W]) for ky
      grid_corr_y = d[1] == 0 ? 1.0 : sin((PI*d[1])/(double)Nmesh)/((PI*d[1])/(double)Nmesh);

      //============================
      // Do the k = 0 mode
      //============================
      d[2] = 0;
      grid_corr_z = 1.0;
      kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
      nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
      if(nk >= 0 && nk < nbins){
        coord = (i*Nmesh+j)*(Nmesh/2+1) + d[2];
        grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fftw_norm_fac;
        pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
        pofk_bin[nk] += 1.0*pofk;
        k_bin[nk]    += 1.0*kmag;
        n_bin[nk]    += 1.0;

        // Mirror around y-axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          coord = (i*Nmesh+(Nmesh-j))*(Nmesh/2+1) + d[2];
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
      grid_corr_z = 2.0 / PI;
      kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
      nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
      if(nk >= 0 && nk < nbins){
        coord = (i*Nmesh+j)*(Nmesh/2+1) + d[2];
        grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fftw_norm_fac;
        pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
        pofk_bin[nk] += 1.0*pofk;
        k_bin[nk]    += 1.0*kmag;
        n_bin[nk]    += 1.0;

        // Mirror around y-axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          coord = (i*Nmesh+(Nmesh-j))*(Nmesh/2+1) + d[2];
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk_bin[nk] += 1.0*pofk;
          k_bin[nk]    += 1.0*kmag;
          n_bin[nk]    += 1.0;
        }
      }

      for (int k = 1; k < (unsigned int) (Nmesh/2); k++) {

        // |kz| component
        d[2] = k;

        // Window-function (Sqrt[W]) for kz
        grid_corr_z = sin((PI*d[2])/(double)Nmesh)/((PI*d[2])/(double)Nmesh);

        //============================
        // Do the general mode
        //============================
        kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
        nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
        if(nk >= 0 && nk < nbins){
          coord = (i*Nmesh+j)*(Nmesh/2+1) + d[2];
          grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fftw_norm_fac;
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk_bin[nk] += 2.0*pofk;
          k_bin[nk]    += 2.0*kmag;
          n_bin[nk]    += 2.0;

          // Mirror around y-axis
          if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
            coord = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + d[2];
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
      pofk_bin[i] = (pofk_bin_all[i] / n_bin_all[i]) * pow( Box, 3);
      if(subtract_shotnoise) pofk_bin[i] -= pow(Box / (double) Nsample, 3);
      k_bin[i]    = (k_bin_all[i]    / n_bin_all[i]) * 2.0 * M_PI / Box;
    }
  }

  // Write output to file
  double znow = 1.0/a - 1.0;
  if(ThisTask == 0){
    
    // Make filename
    char filename[1000];
    int zint = (int) (znow);
    int zfrac = (int)((znow - zint)*1000);
    sprintf(filename, "%s/pofk_%s_z%d.%03d_%s.txt", OutputDir, FileBase, zint, zfrac, label);
    printf("Writing power-spectrum to file: [%s]\n", filename);

    FILE *fp = fopen(filename,"w");
    fprintf(fp,"#  k_bin (h/Mpc)        P(k) (Mpc/h)^3     k_mean_bin (h/Mpc)     Delta = k^3P(k)/2pi^2\n");
    for(int i = 1; i < nbins; i++){
      if(n_bin_all[i] > 0){
        double k_of_bin   = k_from_index(i, kmin, kmax, nbins, bintype);
        double k_mean_bin = k_bin[i];
        double pofk       = pofk_bin[i];
        double Delta      = pofk_bin[i] * k_of_bin * k_of_bin * k_of_bin / 2.0 / M_PI / M_PI;
        fprintf(fp, "%10.5f   %10.5f   %10.5f   %10.5f\n", k_of_bin, pofk, k_mean_bin, Delta);
      }
    }
    fclose(fp);
  }

  // Free memory
  my_free(pofk_bin);
  my_free(pofk_bin_all);
  my_free(n_bin);
  my_free(n_bin_all);
  my_free(k_bin);
  my_free(k_bin_all);
  timer_stop(_PofkComputation);
}

//==========================================================
// Does Cloud-in-Cell assignment where the axis [axis] is 
// transformed to redshift space.
// NB: for scaledependent version of the code we assumeE
// that dDdy[] contains the dD/dy and not DeltaD = Int dD/dy
// over a time-step (as is needed for the time-stepping)
//==========================================================
void PtoMesh_RSD(float_kind *dens, int axis, double A) {
  unsigned int i;
  unsigned int IX, IY, IZ;
  unsigned int IXneigh, IYneigh, IZneigh;
  double X, Y, Z, V;
  double TX, TY, TZ;
  double DX, DY, DZ;
  double scaleBox = (double)Nmesh/Box;
  double WPAR = pow((double)Nmesh / (double)Nsample,3);
  
  // Normalization of velocity 
  double velfac  = (Hubble / A);
  double vnorm   = velfac/(100.0 * A * hubble(A)) * scaleBox;

  if( ! (axis == YAXIS || axis == ZAXIS) ){
    printf("Error: axis [%i] not allowed, putting axis to z = [%i]\n", axis, ZAXIS);
    axis = ZAXIS;
  }

#ifndef SCALEDEPENDENT
  // If not scaledependent growth
  double dDdy  = growth_dDdy(A);
  double dD2dy = growth_dD2dy(A);
#endif

  // Initialize density to -1
  for(i = 0; i < 2 * Total_size; i++) 
    dens[i] = -1.0;
  
  for(i = 0; i < NumPart; i++) {

    // Scale positions to be in [0, Nmesh]
    if(axis == YAXIS){
      X = P[i].Pos[0] * scaleBox;
      Y = P[i].Pos[2] * scaleBox;
      Z = P[i].Pos[1] * scaleBox;
    } else if(axis == ZAXIS){
      X = P[i].Pos[0] * scaleBox;
      Y = P[i].Pos[1] * scaleBox;
      Z = P[i].Pos[2] * scaleBox;
    }
    
    // Transform axis to redshift space
#ifdef SCALEDEPENDENT
    V = UseCOLA == 0 ? P[i].Vel[axis] : P[i].Vel[axis] + (P[i].dDdy[axis] + P[i].dD2dy[axis] );
#else
    V = UseCOLA == 0 ? P[i].Vel[axis] : P[i].Vel[axis] + (P[i].D[axis] * dDdy + P[i].D2[axis] * dD2dy);
#endif

    V *= vnorm;
    Z += V;

    // Periodic BC
    if(Z >= (double) Nmesh) Z -= (double) Nmesh;
    if(Z < 0)               Z += (double) Nmesh;

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
    IX -= Local_x_start;
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
    dens[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]                += TX*TY*TZ;
    dens[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]           += TX*TY*DZ;
    dens[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]           += TX*DY*TZ;
    dens[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]      += TX*DY*DZ;
    dens[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]           += DX*TY*TZ;
    dens[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]      += DX*TY*DZ;
    dens[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]      += DX*DY*TZ;
    dens[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh] += DX*DY*DZ;
  }

  //====================================================================================
  // Copy across the extra slice from the task on the left and add it to the leftmost slice
  // of the task on the right. Skip over tasks without any slices.
  //====================================================================================

  float_kind * temp_density = (float_kind *) my_calloc(2*alloc_slice,sizeof(float_kind));
  ierr = MPI_Sendrecv(&(dens[2*last_slice]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask,0,
      &(temp_density[0]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,LeftTask,0,MPI_COMM_WORLD,&status);
  if (NumPart != 0) {
    for (i = 0; i < 2 * alloc_slice; i++) dens[i] += (temp_density[i] + 1.0);
  }
  my_free(temp_density);

  return;
}
 
//===========================================================================
// Computes the RSD power-spectrum. This requires 2 FFT's for non-scaledependent
// case. For scaledependent (when UseCOLA = 1) we needto have  dD/dy in P.dDdy
// to do this routine. This is the case in the Output(.) routine so it's best 
// to call it there with dDdy_set_in_particles = 1 to save computation. 
// If this is called with dDdy_set_in_particles = 0 when scaledependent then 
// 2 FFT's + temp memory is needed to get the velocity here
//===========================================================================
void compute_RSD_powerspectrum(double A, int dDdy_set_in_particles){
  timer_start(_PofkComputation);

  if(ThisTask == 0)
    printf("Computing the RSD power-spectrum...\n");

#ifdef SCALEDEPENDENT
  float *tmp_buffer_D1 = NULL, *tmp_buffer_D2 = NULL;
  if( UseCOLA && ! dDdy_set_in_particles){
    // Make a copy of the displacment-field [dDdy] which contains [deltaD]
    // and recompute it to give dDdy which we need to set velocites below
    tmp_buffer_D1 = my_malloc(sizeof(float) * NumPart * 3);
    tmp_buffer_D2 = my_malloc(sizeof(float) * NumPart * 3);
    for(int i = 0; i < NumPart; i++){
      for(int axes = 0; axes < 3; axes++){
        tmp_buffer_D1[3*i + axes] = P[i].dDdy[axes];
        tmp_buffer_D2[3*i + axes] = P[i].dD2dy[axes];
      }
    }

    // Compute dDdy
    if(ThisTask == 0) printf("We need to compute velocities for particles (2 FFT's)...\n");
    double As = aexp_global;
    assign_displacment_field_to_particles(As, As, As, FIELD_dDdy, LPT_ORDER_ONE);
    assign_displacment_field_to_particles(As, As, As, FIELD_dDdy, LPT_ORDER_TWO);
  }
#endif

  // Allocate memory for density grid 
  float_kind *dens = my_malloc(2 * Total_size * sizeof(float_kind));
  complex_kind *dens_k = (complex_kind *) dens;
  plan_kind densplan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh, dens, dens_k, MPI_COMM_WORLD, FFTW_ESTIMATE);

  // Do y-axis
  if(ThisTask == 0) printf("FFT redshift space density (LOS = y-axis)...\n");
  PtoMesh_RSD(dens, YAXIS, A);
  my_fftw_execute(densplan);
  struct RSD_Pofk_Data pofk_data_y;
  bin_up_RSD_power_spectrum(dens_k, &pofk_data_y);

  // Do z-axis
  if(ThisTask == 0) printf("FFT redshift space density (LOS = z-axis)...\n");
  PtoMesh_RSD(dens, ZAXIS, A);
  my_fftw_execute(densplan);
  struct RSD_Pofk_Data pofk_data_z;
  bin_up_RSD_power_spectrum(dens_k, &pofk_data_z);

  // Doing the x-axis requires communication as we will need to move particles to the left task
  // so we skip this for now...

  // Make average of the two axes and output
  if(ThisTask == 0){

    // Make filename
    char filename[1000];
    double znow = 1.0/A - 1.0;
    int zint = (int) (znow);
    int zfrac = (int)((znow - zint)*1000);
    sprintf(filename, "%s/pofk_RSD_%s_z%d.%03d.txt", OutputDir, FileBase, zint, zfrac);
    printf("Writing RSD power-spectrum to file: [%s]\n", filename);

    FILE *fp = fopen(filename,"w");
    fprintf(fp,"#  k (h/Mpc)      P0 (Mpc/h)^3      P2 (Mpc/h)^3      P4 (Mpc/h)^3     sigma0     sigma2     sigma4\n");
    for(int i = 0; i < pofk_data_z.nbins; i++){
      if(pofk_data_y.n[i] > 0 && pofk_data_y.k[i] > 0.0){
        double k    = pofk_data_y.k[i];
        double P0   = (pofk_data_y.P0[i] + pofk_data_z.P0[i])/2.0;
        double P2   = (pofk_data_y.P2[i] + pofk_data_z.P2[i])/2.0;
        double P4   = (pofk_data_y.P4[i] + pofk_data_z.P4[i])/2.0;
        double err0 = fabs(pofk_data_y.P0[i] - pofk_data_z.P0[i]) / sqrt(2.0);
        double err2 = fabs(pofk_data_y.P2[i] - pofk_data_z.P2[i]) / sqrt(2.0);
        double err4 = fabs(pofk_data_y.P4[i] - pofk_data_z.P4[i]) / sqrt(2.0);
        fprintf(fp, "%10.5f   %10.5f   %10.5f   %10.5f   %10.5f   %10.5f   %10.5f\n", k, P0, P2, P4, err0, err2, err4);
      }
    }
    fclose(fp);
  }

  // Clean up
  my_free(pofk_data_y.n );
  my_free(pofk_data_y.k );
  my_free(pofk_data_y.P0);
  my_free(pofk_data_y.P2);
  my_free(pofk_data_y.P4);

  my_free(pofk_data_z.n );
  my_free(pofk_data_z.k );
  my_free(pofk_data_z.P0);
  my_free(pofk_data_z.P2);
  my_free(pofk_data_z.P4);
  
  my_free(dens);
  my_fftw_destroy_plan(densplan);

#ifdef SCALEDEPENDENT
  if( UseCOLA && ! dDdy_set_in_particles){
    // Copy back copy of [deltaD] to [dDdy] and free up memory
    for(int i = 0; i < NumPart; i++){
      for(int axes = 0; axes < 3; axes++){
        P[i].dDdy[axes]  = tmp_buffer_D1[3*i + axes];
        P[i].dD2dy[axes] = tmp_buffer_D2[3*i + axes];
      }
    }
    my_free(tmp_buffer_D1);
    my_free(tmp_buffer_D2);
  }
#endif

  timer_stop(_PofkComputation);
}

//===========================================================================
// Compute P0, P2, P4 RSD power-spectra
// Assumes dens_k has the fourier transform of the density field
//===========================================================================
void bin_up_RSD_power_spectrum(complex_kind *dens_k, struct RSD_Pofk_Data *pofk_data){
  unsigned int coord;
  double pofk, mu2, kmag, grid_corr_x, grid_corr_y, grid_corr_z, grid_corr = 1.0;
  int d[3], nk;

  //=======================================================
  // Define the binning. Bintype: 0 = Linear; 1 = Logspaced
  // Every integer k-mode: kmin = 0, kmax = Nmesh/2, 
  // nbins = Nmesh/2+1 and bintype = 0
  //=======================================================

#ifdef COMPUTE_POFK

  // Get parameters from parameterfile
  int nbins              = pofk_nbins;                     // Number of bins in k
  int bintype            = pofk_bintype;;                  // Linear [0] Logarithmic [1]
  int subtract_shotnoise = pofk_subtract_shotnoise;        // Shot-noise subtraction
  double kmin            = pofk_kmin * Box / (2.0 * M_PI); // Min 'integer' wave-number
  double kmax            = pofk_kmax * Box / (2.0 * M_PI); // Max 'integer' wave-number

#else

  // Fiducial parameters
  int nbins              = Nmesh/2;
  int bintype            = LINEAR_SPACING;
  int subtract_shotnoise = 1;
  double kmin            = 0.0;
  double kmax            = Nmesh/2.0;

#endif
  
  // Sanity checks. Adjust parameters to sane values
  adjust_pofk_parameters(&nbins, &bintype, &subtract_shotnoise, &kmin, &kmax);

  double *pofk0_bin     = my_malloc(sizeof(double) * nbins);
  double *pofk2_bin     = my_malloc(sizeof(double) * nbins);
  double *pofk4_bin     = my_malloc(sizeof(double) * nbins);
  double *pofk0_bin_all = my_malloc(sizeof(double) * nbins);
  double *pofk2_bin_all = my_malloc(sizeof(double) * nbins);
  double *pofk4_bin_all = my_malloc(sizeof(double) * nbins);
  double *n_bin         = my_malloc(sizeof(double) * nbins);
  double *n_bin_all     = my_malloc(sizeof(double) * nbins);
  double *k_bin         = my_malloc(sizeof(double) * nbins);
  double *k_bin_all     = my_malloc(sizeof(double) * nbins);

  for(int i = 0; i < nbins; i++){
    pofk0_bin[i] = pofk0_bin_all[i] = 0.0;
    pofk2_bin[i] = pofk2_bin_all[i] = 0.0;
    pofk4_bin[i] = pofk4_bin_all[i] = 0.0;
    n_bin[i]     = n_bin_all[i]     = 0.0;
    k_bin[i]     = k_bin_all[i]     = 0.0;
  }
  
  // FFT normalization factor for |density(k)|^2
  double fftw_norm_fac = 1.0 / pow( (double) Nmesh, 6);
        
  // Loop over all modes and add up RSD multipoles P0, P2, P4
  for (int i = 0; i < Local_nx; i++) {
    int iglobal = i + Local_x_start;

    // |kx| component
    d[0] = iglobal > Nmesh/2 ? Nmesh - iglobal : iglobal;

    // Window-function (Sqrt[W]) for kx
    grid_corr_x = d[0] == 0 ? 1.0 : sin((PI*d[0])/(double)Nmesh)/((PI*d[0])/(double)Nmesh);

    for (int j = 0 ; j < (unsigned int)(Nmesh/2+1); j++) {

      // |ky| component
      d[1] = j;

      // Window-function (Sqrt[W]) for ky
      grid_corr_y = d[1] == 0 ? 1.0 : sin((PI*d[1])/(double)Nmesh)/((PI*d[1])/(double)Nmesh);

      //============================
      // Do the k = 0 mode
      //============================
      d[2] = 0;
      grid_corr_z = 1.0;
      kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
      mu2 = d[2]*d[2]/(kmag*kmag);
      if(iglobal == 0 && j == 0) mu2 = 0.0;
      nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
      if(nk >= 0 && nk < nbins){
        coord = (i*Nmesh+j)*(Nmesh/2+1) + d[2];
        grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fftw_norm_fac;
        pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
        pofk0_bin[nk] += 1.0*pofk;
        pofk2_bin[nk] += 1.0*pofk*mu2;
        pofk4_bin[nk] += 1.0*pofk*mu2*mu2;
        k_bin[nk]     += 1.0*kmag;
        n_bin[nk]     += 1.0;

        // Mirror around y-axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          coord = (i*Nmesh+(Nmesh-j))*(Nmesh/2+1) + d[2];
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk0_bin[nk] += 1.0*pofk;
          pofk2_bin[nk] += 1.0*pofk*mu2;
          pofk4_bin[nk] += 1.0*pofk*mu2*mu2;
          k_bin[nk]     += 1.0*kmag;
          n_bin[nk]     += 1.0;
        }
      }

      //============================
      // Do the k = N/2 mode
      //============================
      d[2] = Nmesh/2;
      grid_corr_z = 2.0 / PI;
      kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
      mu2 = (d[2]*d[2]/kmag/kmag);
      nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
      if(nk >= 0 && nk < nbins){
        coord = (i*Nmesh+j)*(Nmesh/2+1) + d[2];
        grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fftw_norm_fac;
        pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
        pofk0_bin[nk] += 1.0*pofk;
        pofk2_bin[nk] += 1.0*pofk*mu2;
        pofk4_bin[nk] += 1.0*pofk*mu2*mu2;
        k_bin[nk]     += 1.0*kmag;
        n_bin[nk]     += 1.0;

        // Mirror around y-axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          coord = (i*Nmesh+(Nmesh-j))*(Nmesh/2+1) + d[2];
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk0_bin[nk] += 1.0*pofk;
          pofk2_bin[nk] += 1.0*pofk*mu2;
          pofk4_bin[nk] += 1.0*pofk*mu2*mu2;
          k_bin[nk]     += 1.0*kmag;
          n_bin[nk]     += 1.0;
        }
      }

      for (int k = 1; k < (unsigned int) (Nmesh/2); k++) {

        // |kz| component
        d[2] = k;

        // Window-function (Sqrt[W]) for kz
        grid_corr_z = sin((PI*d[2])/(double)Nmesh)/((PI*d[2])/(double)Nmesh);

        //============================
        // Do the general mode
        //============================
        kmag = sqrt( d[0] * d[0] + d[1] * d[1] + d[2] * d[2] );
        mu2 = (d[2]*d[2]/kmag/kmag);
        nk = pofk_bin_index(kmag, kmin, kmax, nbins, bintype);
        if(nk >= 0 && nk < nbins){
          coord = (i*Nmesh+j)*(Nmesh/2+1) + d[2];
          grid_corr = 1.0 / pow(grid_corr_x * grid_corr_y * grid_corr_z, 4.0) * fftw_norm_fac;
          pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
          pofk0_bin[nk] += 2.0*pofk;
          pofk2_bin[nk] += 2.0*pofk*mu2;
          pofk4_bin[nk] += 2.0*pofk*mu2*mu2;
          k_bin[nk]     += 2.0*kmag;
          n_bin[nk]     += 2.0;

          // Mirror around y-axis
          if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
            coord = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + d[2];
            pofk = (dens_k[coord][0] * dens_k[coord][0] + dens_k[coord][1] * dens_k[coord][1]) * grid_corr;
            pofk0_bin[nk] += 2.0*pofk;
            pofk2_bin[nk] += 2.0*pofk*mu2;
            pofk4_bin[nk] += 2.0*pofk*mu2*mu2;
            k_bin[nk]     += 2.0*kmag;
            n_bin[nk]     += 2.0;
          }
        }
      }
    }
  }

  // Communicate
  MPI_Allreduce(pofk0_bin, pofk0_bin_all, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pofk2_bin, pofk2_bin_all, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pofk4_bin, pofk4_bin_all, nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(n_bin,     n_bin_all,     nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(k_bin,     k_bin_all,     nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Normalize and subtract shotnoise
  for(int i = 0; i < nbins; i++){
    
    if(n_bin_all[i] > 0){
      pofk0_bin[i] = (pofk0_bin_all[i] / n_bin_all[i]) * pow( Box, 3);
      pofk2_bin[i] = (pofk2_bin_all[i] / n_bin_all[i]) * pow( Box, 3);
      pofk4_bin[i] = (pofk4_bin_all[i] / n_bin_all[i]) * pow( Box, 3);

      // From powers of mu to Legendre polynomial combinations (NB: must be done in correct order)
      pofk4_bin[i] = 9.0*(35.0*pofk4_bin[i] - 30.0*pofk2_bin[i] + 3.0*pofk0_bin[i])/8.0;
      pofk2_bin[i] = 5.0*(3.0*pofk2_bin[i] - 1.0*pofk0_bin[i])/2.0;
      pofk0_bin[i] = 1.0*pofk0_bin[i]; 

      if(subtract_shotnoise){
        pofk0_bin[i] -= pow(Box / (double) Nsample, 3);
        pofk2_bin[i] -= pow(Box / (double) Nsample, 3);
        pofk4_bin[i] -= pow(Box / (double) Nsample, 3);
      }

      k_bin[i] = k_from_index(i, kmin, kmax, nbins, bintype);

      // Alternative (use mean of all k-modes we added to the bin):
      //k_bin[i]    = (k_bin_all[i]    / n_bin_all[i]) * 2.0 * M_PI / Box;
    }

    n_bin[i] = n_bin_all[i];
  }

  // Copy over the results
  pofk_data->nbins = nbins;
  pofk_data->n  = my_malloc(sizeof(double)*nbins);
  pofk_data->k  = my_malloc(sizeof(double)*nbins);
  pofk_data->P0 = my_malloc(sizeof(double)*nbins);
  pofk_data->P2 = my_malloc(sizeof(double)*nbins);
  pofk_data->P4 = my_malloc(sizeof(double)*nbins);
  for(int i = 0; i < nbins; i++){
    pofk_data->n[i]  = n_bin[i];
    pofk_data->k[i]  = k_bin[i];
    pofk_data->P0[i] = pofk0_bin[i];
    pofk_data->P2[i] = pofk2_bin[i];
    pofk_data->P4[i] = pofk4_bin[i];
  }

  // Free memory
  my_free(pofk0_bin);
  my_free(pofk2_bin);
  my_free(pofk4_bin);
  my_free(pofk0_bin_all);
  my_free(pofk2_bin_all);
  my_free(pofk4_bin_all);
  my_free(n_bin);
  my_free(n_bin_all);
  my_free(k_bin);
  my_free(k_bin_all);
}

//===========================================================
// Check and adjust P(k) settings if they don't make sense
//===========================================================
void adjust_pofk_parameters(int *nbins, int *bintype, int *subtract_shotnoise, double *kmin, double *kmax){
  
  // Sanity checks. Adjust parameters to sane values
  if( ! (*bintype == LINEAR_SPACING || *bintype == LOG_SPACING) ){
    if(ThisTask == 0) printf("pofk: Unknown bintype [%i]. Will use linear-spacing\n", *bintype);
    *bintype = LINEAR_SPACING;
  }
  if(! (*subtract_shotnoise == 0 || *subtract_shotnoise == 1)){
    if(ThisTask == 0) printf("pofk: Unknown shotnoise option [%i] != 0 or 1. Setting it to [1] = true\n", *subtract_shotnoise);
    *subtract_shotnoise = 1;
  }
  if( *nbins <= 0 ){
    if(ThisTask == 0) printf("pofk: Nbins = [%i] < 0. Using nbins = [%i] instead\n", *nbins, Nmesh);
    *nbins = Nmesh;
  }
  int adjustrange = 0;
  if(*kmax <= *kmin){
    if(*bintype == LINEAR_SPACING) *kmin = 0.0;
    if(*bintype == LOG_SPACING)    *kmin = 1.0;
    *kmax = (double) Nmesh;
    adjustrange = 1;
  }
  if(*kmin < 0.0){
    if(*bintype == LINEAR_SPACING) *kmin = 0.0;
    if(*bintype == LOG_SPACING)    *kmin = 1.0;
    adjustrange = 1;
  }
  if( (*bintype == LOG_SPACING) && (*kmin == 0.0) ){
    *kmin = 1.0;
    adjustrange = 1;
  }
  if(*kmax > sqrt(3.0) * (double) Nmesh){
    *kmax = (double) Nmesh; 
    adjustrange = 1;
  }
  if(ThisTask == 0 && adjustrange) 
    printf("pofk: We had to adjust the k-range kmin = [%5.3f] kmax = [%5.3f]  h/Mpc \n", 2.0 * M_PI / Box * (*kmin), 2.0 * M_PI / Box * (*kmax));

#ifdef COMPUTE_POFK
  // Update the values so we don't have to adjust it again
  pofk_kmin    = (*kmin) * 2.0 * M_PI / Box;
  pofk_kmax    = (*kmax) * 2.0 * M_PI / Box;
  pofk_nbins   = *nbins;
  pofk_bintype = *bintype;
  pofk_subtract_shotnoise = *subtract_shotnoise;
#endif

}

