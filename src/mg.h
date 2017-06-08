#ifndef INCLUDEMG
#define INCLUDEMG

//==========================================================================//
//                                                                          //
//  MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017           //
//                                                                          //
// This file contains all the modified gravity modifications to the         //
// MG-PICOLA code needed to solve the scalar field equation including       //
// screening                                                                //
//                                                                          //
//==========================================================================//

//=============================================
// Takes the k-space complex grid densityk and 
// divides by Laplacian and stores it in phinewton(k)
//
// Needed for chameleon type models
//=============================================

void DivideByLaplacian(complex_kind *densityk, complex_kind *phinewtonk){
  int iglobal, kmin;
  double RK, KK;

  // Normalization of the density as when *= -1/k^2 and FFTd gives us Phi(x) in correct units
  // The Phi we compute here is the one that satisfy D_x^2 Phi = 4 pi G a^2 rho(a) delta = 1.5 Omega/a H0^2 delta
  double normfactor = 1.0/pow((double)Nmesh,3);
  normfactor *= 1.5 * Omega / aexp_global * pow (Box / INVERSE_H0_MPCH / (2.0 * PI) , 2);

  // We need global values for i as opposed to local values
  // Same goes for anything that relies on i (such as RK). 
  for (int i = 0; i < Local_nx; i++) {
    iglobal = i + Local_x_start;
    for (int j = 0; j < (unsigned int)(Nmesh/2+1); j++) {
      kmin = 0;
      if ((iglobal == 0) && (j == 0)) {
        phinewtonk[0][0] = phinewtonk[0][1] = 0.0;
        kmin = 1;
      }
      for (int k = kmin; k < (unsigned int)(Nmesh/2+1); k++) {
        unsigned int ind = (i*Nmesh + j)*(Nmesh/2+1) + k;
        if (iglobal > Nmesh/2) {
          RK    = (double)(k*k + (Nmesh-iglobal)*(Nmesh-iglobal) + j*j);
        } else {
          RK    = (double)(k*k + iglobal*iglobal + j*j);
        }
        KK = -1.0/RK;

        // Divide by laplacian
        phinewtonk[ind][0] = (normfactor*densityk[ind][0]*KK);
        phinewtonk[ind][1] = (normfactor*densityk[ind][1]*KK);

        // Do the mirror along the y axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          ind = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + k;

          // Divide by laplacian
          phinewtonk[ind][0] = (normfactor*densityk[ind][0]*KK);
          phinewtonk[ind][1] = (normfactor*densityk[ind][1]*KK);
        }
      }
    }
  }
}

//=============================================
// Change density_eff(k) to phi(k) which is 
// normalized such that it can just be added
// to the Newtonian force
//
// Needed for chameleon type models
//=============================================

void EffDensitykToPhiofk(complex_kind* P3D_densityeffk, complex_kind *P3D_phik){
  int iglobal, kmin;
  double RK, KK;

  // Normalization factors
  double normfactor = coupling_function(aexp_global);
  double massterm2 = pow(aexp_global, 2) * mass2_of_a(aexp_global) / pow( (2.0 * PI) * INVERSE_H0_MPCH / Box, 2);

  // We need global values for i as opposed to local values
  // Same goes for anything that relies on i (such as RK). 
  for (int i = 0; i < Local_nx; i++) {
    iglobal = i + Local_x_start;
    for (int j = 0; j < (unsigned int)(Nmesh/2+1); j++) {
      kmin = 0;
      if ((iglobal == 0) && (j == 0)) {
        P3D_phik[0][0] = P3D_phik[0][1] = 0.0;
        kmin = 1;
      }
      for (int k = kmin; k < (unsigned int)(Nmesh/2+1); k++) {
        unsigned int ind = (i*Nmesh + j)*(Nmesh/2+1) + k;
        if (iglobal > Nmesh/2) {
          RK    = (double)(k*k + (Nmesh-iglobal)*(Nmesh-iglobal) + j*j);
        } else {
          RK    = (double)(k*k + iglobal*iglobal + j*j) ;
        }
        KK = RK / (RK + massterm2);

        // Normalize
        P3D_phik[ind][0] = (normfactor * P3D_densityeffk[ind][0]*KK);
        P3D_phik[ind][1] = (normfactor * P3D_densityeffk[ind][1]*KK);

        // Do the mirror along the y axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          ind = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + k;

          // Normalize
          P3D_phik[ind][0] = (normfactor * P3D_densityeffk[ind][0]*KK);
          P3D_phik[ind][1] = (normfactor * P3D_densityeffk[ind][1]*KK);
        }
      }
    }
  }
}

//=====================================================================
// The main driver for any MG models that only has a modified Geff(a)
// so that we don't need to compute any fifth-force potential, we 
// simply rescale the density array (the source for the force)
//
// Assuming we have density(k) in P3D when starting this routine
// as CopyDensityArray has been called in PtoMesh
//=====================================================================

void ComputeFifthForce_TimeDepGeffModels(){
  if(ThisTask == 0)
    printf("\n===> Multiplying with Geff(a) \n");
  
  double Geffective = GeffoverG(aexp_global, 0.0);

  for(unsigned int j  = 0; j < Total_size; j++){
    P3D[j][0] *= Geffective;
    P3D[j][1] *= Geffective;
  }
}

//=====================================================================
// The main driver for any scalar tensor gravity defined by m(a) and
// beta(a) like f(R) gravity
//
// Assuming we have density(k) in P3D when starting this routine
// as CopyDensityArray has been called in PtoMesh
//=====================================================================

void ComputeFifthForce_PotentialScreening(){
  if(ThisTask == 0)
    printf("\n===> Computing modified gravity potential\n");

  // When no screening we simply just need to call EffDensitykToPhiofk
  if( !include_screening ){
    my_fftw_execute(plan_mg_phik);
    EffDensitykToPhiofk(P3D, P3D_mgarray_two);
    return;
  }

  //=====================================================================
  // Loop over k-grid and divide by laplacian
  //=====================================================================
  DivideByLaplacian(P3D, P3D_mgarray_one);

  //=====================================================================
  // Fourier-transform P3D_mgarray_one to get Phi(x).
  // After this is done we have Phi(x) in mgarray_one
  //=====================================================================
  if(ThisTask == 0)
    printf("===> Fourier transforming to get Phi(x)\n");
  my_fftw_execute(plan_mg_phinewton);

  //=====================================================================
  // Compute density_effective(x) and store it in mgarray_two
  //=====================================================================
  for(unsigned int j  = 0; j < 2*Total_size; j++){
    mgarray_two[j] *= screening_factor_potential(aexp_global, mgarray_one[j]);
  }

  //=====================================================================
  // Fourier-transform to get density_eff(k) stored in P3D_mgarray_two
  //=====================================================================
  if(ThisTask == 0)
    printf("===> Fourier transforming to get effective density(k)\n");
  my_fftw_execute(plan_mg_phik);

  //=====================================================================
  // Compute force potential phi(k) = density_eff(k) * k^2 / (k^2 + m^2)
  //=====================================================================
  EffDensitykToPhiofk(P3D_mgarray_two, P3D_mgarray_two);
}

//=====================================================================
// The main driver for DGP like models that screen depending on density
// Assuming we have density(k) in P3D when starting this routine
// as CopyDensityArray has been called in PtoMesh
//=====================================================================

void ComputeFifthForce_DensityScreening(){
  double Rsmooth  = 1.0; // Smoothing radius in Mpc/h
  double coupling = coupling_function(aexp_global);

#ifdef DGPGRAVITY

  // Use the smoothing radius set in the parameterfile
  Rsmooth = Rsmooth_global; 

#endif

  //=====================================================================
  // If include_screening is not active just assign phik to be 
  // the density(k) * coupling
  //=====================================================================
  if( ! include_screening ){

    for(unsigned int i = 0; i < 2 * Total_size; i++){
      mgarray_two[i] = density[i] * coupling;
    }
    return;
  }

  //=====================================================================
  // For testing: Check that the copy we took of the density field works
  // if(ThisTask == 0) printf("     Density field before smoothing:\n");
  // check_real_space_grid(mgarray_two);
  //=====================================================================

  //=====================================================================
  // Compute smoothed density by multiplying by filter in k-space
  //=====================================================================
  SmoothDensityField(P3D, P3D_mgarray_one, Rsmooth);

  //=====================================================================
  // Transform to real-space. After this we should have the 
  // smoothed density field in mgarray_one
  //=====================================================================
  my_fftw_execute(plan_mg_phinewton);

  //=====================================================================
  // Compute density_effective(x) and store it in mgarray_two
  //=====================================================================
  double maxscreen = 0.0, minscreen = 1e100, avgscreen = 0.0;
  for(unsigned int j = 0; j < 2 * Total_size; j++){
    double screenfac = screening_factor_density(aexp_global, mgarray_one[j]);
    mgarray_two[j] *= coupling * screenfac;
    if(screenfac > maxscreen) maxscreen = screenfac;
    if(screenfac < minscreen) minscreen = screenfac;
    avgscreen += screenfac;

    // Alternative: use smoothed density field also for force, not just for the screening factor
    // mgarray_two[j] = coupling * mgarray_one[j] * screening_factor_density(aexp_global, mgarray_one[j]);
  }
  avgscreen /= (double)(2 * Total_size);
  if(ThisTask == 0){
    printf("===> For CPU[0] we have AvgScreenFac = %8.3f  MaxScreenFac = %8.3f  MinScreenFac = %8.3f\n", avgscreen, maxscreen, minscreen);
  }

  //=====================================================================
  // Transform effective density to k-space
  //=====================================================================
  my_fftw_execute(plan_mg_phik);
}

//=============================================
// Smooth the density field by multiplying
// the Fourier transform of the density
// by a k-space smoothing filter of size Rsmooth
//=============================================

void SmoothDensityField(complex_kind *densityk, complex_kind *densityk_smooth, double Rsmooth){
  int iglobal;
  double RK;

  // Normalization of the density as when *= -1/k^2 and FFTd gives us Phi(x) in correct units
  double normfactor = 1.0/pow((double)Nmesh,3);

  if(ThisTask == 0) {
    printf("\n===> Compute smoothing factor for wavenumbers in k_min = %f h/Mpc  ->  k_max = %f h/Mpc\n",  2.0 * PI / Box, 2.0 * PI / Box * Nmesh * sqrt(3.0));
    printf(  "===> The filter supression factor is f(kmin) = %f  and  f(kmax) = %f\n",  smoothing_filter(2.0 * PI/ Box * Rsmooth),  smoothing_filter(2.0 * PI/ Box * Rsmooth * Nmesh / 2.0) );
  }

  // We need global values for i as opposed to local values
  // Same goes for anything that relies on i (such as RK). 
  for (int i = 0; i < Local_nx; i++) {
    iglobal = i + Local_x_start;
    for (int j = 0; j < (unsigned int)(Nmesh/2+1); j++) {
      for (int k = 0; k < (unsigned int)(Nmesh/2+1); k++) {
        unsigned int ind = (i*Nmesh + j)*(Nmesh/2+1) + k;
        if (iglobal > Nmesh/2) {
          RK    = (double)(k*k + (Nmesh-iglobal)*(Nmesh-iglobal) + j*j);
        } else {
          RK    = (double)(k*k + iglobal*iglobal + j*j) ;
        }

        // Smooth the density field
        double smooth = smoothing_filter( sqrt(RK) * 2.0 * PI/ Box * Rsmooth) * normfactor;
        densityk_smooth[ind][0] = densityk[ind][0] * smooth;
        densityk_smooth[ind][1] = densityk[ind][1] * smooth;

        // Do the mirror along the y axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          ind = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + k;

          // Smooth the density field
          densityk_smooth[ind][0] = densityk[ind][0] * smooth;
          densityk_smooth[ind][1] = densityk[ind][1] * smooth;
        }
      }
    }
  }
}

//=============================================
// Smoothing filter for the density field
// kR is wavenumber times smooth radius
// Needed for models that have a screening factor
// depending on density
//=============================================

double smoothing_filter(double kR){
#if defined(GAUSSIANFILTER)

  return exp(-0.5 * kR * kR);

#elif defined(TOPHATFILTER)

  if(kR < 1e-5) return 1.0; 
  return 3.0/(kR*kR*kR) * (sin(kR) - kR * cos(kR));

#elif defined(SHARPKFILTER)

  if(kR < 1.0) return 1.0;
  return 0.0;

#else

  return 1.0;

#endif
}

//=============================================
// This method computed the mean, rms, max 
// and min of a real-space grid
//=============================================

void check_real_space_grid(float_kind *grid){
  double maxval = -1e100;
  double minval = 1e100;
  double avg    = 0.0;
  double rms    = 0.0;
  double ncells = pow((double) Nmesh,3);;

  for(int ix = 0; ix < Local_nx; ix++){
    for(int iy = 0; iy < Nmesh; iy++){
      for(int iz = 0; iz < Nmesh; iz++){
        unsigned int i = (ix * Nmesh + iy)*2*(Nmesh/2 + 1) + iz;
        avg += grid[i];
        rms += grid[i]*grid[i];
        if(grid[i] < minval) minval = grid[i];
        if(grid[i] > maxval) maxval = grid[i];
      }
    }
  }

  // Communicate
  ierr = MPI_Allreduce(MPI_IN_PLACE, &minval, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE, &maxval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE, &avg,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ierr = MPI_Allreduce(MPI_IN_PLACE, &rms,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  avg /= ncells;
  rms  = sqrt(rms / ncells);

  if(ThisTask == 0) printf("     Checking grid: minval = %f   maxval  = %f  avg = %f  rms = %f\n", minval, maxval, avg, rms);
}

//=============================================
// Take a copy of density grid. 
// We need this to compute effective density
//=============================================

void CopyDensityArray(){
  if(allocate_mg_arrays){
    for(unsigned int i = 0; i < 2*Total_size; i++)
      mgarray_two[i] = density[i];
  }
}

//=============================================
// Allocate all the arrays and the
// FFT plans we need
//=============================================

void AllocateMGArrays(){
  if(allocate_mg_arrays){
    mgarray_one       = (float_kind *)   my_malloc(2 * Total_size * sizeof(float_kind));
    mgarray_two       = (float_kind *)   my_malloc(2 * Total_size * sizeof(float_kind));
    P3D_mgarray_one   = (complex_kind *) mgarray_one;
    P3D_mgarray_two   = (complex_kind *) mgarray_two;
    plan_mg_phinewton = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, P3D_mgarray_one, mgarray_one,     MPI_COMM_WORLD, FFTW_ESTIMATE);
    plan_mg_phik      = my_fftw_mpi_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh, mgarray_two,     P3D_mgarray_two, MPI_COMM_WORLD, FFTW_ESTIMATE);
  }
}

//=============================================
// Free up all the arrays we have allocated
// and delete FFT plans we have made
//=============================================

void FreeMGArrays(){
  if(allocate_mg_arrays){
    my_free(mgarray_one);
    my_free(mgarray_two);
    my_fftw_destroy_plan(plan_mg_phinewton);
    my_fftw_destroy_plan(plan_mg_phik);
  }
}

//==========================================================
// For models where density depends on |DPhi|^2
// Just added some routines for this for completeness
// Not tested and there are faster ways to compute this
// (compute gradient in real-space is probably better)
//==========================================================

void ComputeFifthForce_GradientScreening(){

  // Make temp array to store density
  float_kind *density_temp = my_malloc(2 * sizeof(float_kind) * Total_size);
  for(unsigned int i = 0; i < 2*Total_size; i++) density_temp[i] = mgarray_two[i];

  for(unsigned int i = 0; i < 2*Total_size; i++) mgarray_two[i] = 0.0;
  
  Density_to_DPhiNewtonk(P3D, P3D_mgarray_one, 0);
  my_fftw_execute(plan_mg_phinewton);
  for(unsigned int i = 0; i < 2*Total_size; i++) mgarray_two[i] += mgarray_one[i]*mgarray_one[i];
  
  Density_to_DPhiNewtonk(P3D, P3D_mgarray_one, 1);
  my_fftw_execute(plan_mg_phinewton);
  for(unsigned int i = 0; i < 2*Total_size; i++) mgarray_two[i] += mgarray_one[i]*mgarray_one[i];
  
  Density_to_DPhiNewtonk(P3D, P3D_mgarray_one, 2);
  my_fftw_execute(plan_mg_phinewton);
  for(unsigned int i = 0; i < 2*Total_size; i++) mgarray_two[i] += mgarray_one[i]*mgarray_one[i];

  // We now have (DPhi)^2 in units of (h/Mpc)^2 in mgarray_two so we can compute
  // effective density
  for(unsigned int i = 0; i < 2*Total_size; i++) 
    mgarray_two[i] = coupling_function(aexp_global) * density_temp[i] * screening_factor_gradient(aexp_global, mgarray_one[i]);
  my_free(density_temp);

  my_fftw_execute(plan_mg_phik);
}

//=============================================================================
// For models where screening depend on |DPhi|^2 we here transform from delta(k)
// to [ D Phi ]_axes(k) which we fourier transformed is in units of 1/Boxsize
//=============================================================================
void Density_to_DPhiNewtonk(complex_kind *densityk, complex_kind *DPhi_i, int axes){
  int iglobal, kmin;
  double RK, KK;

  // Normalization of the density as when *= -1/k^2 and FFTd gives us Phi(x) in correct units
  // The Phi we compute here is the one that satisfy D_x^2 Phi = 4 pi G a^2 rho(a) delta = 1.5 Omega/a H0^2 delta
  double normfactor = 1.0/pow((double)Nmesh,3);
  normfactor *= 1.5 * Omega / aexp_global * pow (Box / INVERSE_H0_MPCH / (2.0 * PI) , 2) * 2.0 * M_PI / Box;

  // We need global values for i as opposed to local values
  // Same goes for anything that relies on i (such as RK). 
  for (int i = 0; i < Local_nx; i++) {
    iglobal = i + Local_x_start;
    for (int j = 0; j < (unsigned int)(Nmesh/2+1); j++) {
      kmin = 0;
      if ((iglobal == 0) && (j == 0)) {
        DPhi_i[0][0] = DPhi_i[0][1] = 0.0;
        kmin = 1;
      }
      for (int k = kmin; k < (unsigned int)(Nmesh/2+1); k++) {
        unsigned int ind = (i*Nmesh + j)*(Nmesh/2+1) + k;
        double d[3];
        if (iglobal > Nmesh/2) {
          RK    = (double)(k*k + (Nmesh-iglobal)*(Nmesh-iglobal) + j*j);
          d[0] = Nmesh - iglobal;
          d[1] = j;
          d[2] = k;
        } else {
          RK    = (double)(k*k + iglobal*iglobal + j*j);
          d[0] = iglobal;
          d[1] = j;
          d[2] = k;
        }
        KK = -1.0/RK;

        // Divide by laplacian and multiply by i k_axes
        DPhi_i[ind][0] = -(normfactor*densityk[ind][1]*KK) * d[axes];
        DPhi_i[ind][1] =  (normfactor*densityk[ind][0]*KK) * d[axes];

        // Do the mirror along the y axis
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          ind = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + k;

          d[1] = Nmesh - j;

          // Divide by laplacian and multiply by i k_axes
          DPhi_i[ind][0] = -(normfactor*densityk[ind][1]*KK * d[axes]);
          DPhi_i[ind][1] =  (normfactor*densityk[ind][0]*KK * d[axes]);
        }
      }
    }
  }
}

#endif
