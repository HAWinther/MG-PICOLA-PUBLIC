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

//=========================================================================================//
//This file contains some additional routines (parallel and serial) needed for any PM code //
//=========================================================================================//

#include "vars.h"
#include "proto.h"
#include "msg.h"
#include "timer.h"
#include "mg.h"            

//=================================================================
// A master routine called from main.c to calculate the acceleration
//=================================================================
void GetDisplacements(void) {

  //=================================================================
  // First we check whether all the particles are on the correct processor after the last time step/
  // original 2LPT displacement and move them if not
  //=================================================================
  if (ThisTask == 0) printf("Moving particles across task boundaries...\n");
  MoveParticles();

#ifdef MEMORY_MODE
  density = my_malloc(2*Total_size*sizeof(float_kind));
  P3D     = (complex_kind *) density;
  plan    = my_fftw_mpi_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh, density, P3D, MPI_COMM_WORLD, FFTW_ESTIMATE);
  if(modified_gravity_active) AllocateMGArrays();
#endif

  //=================================================================
  // Then we do the Cloud-in-Cell assignment to get the density grid and FFT it.  
  //=================================================================
  if (ThisTask == 0) printf("Calculating density using Cloud-in-Cell...\n");
  PtoMesh();

  if(modified_gravity_active) ComputeFifthForce();

#ifdef MEMORY_MODE
  N11  = my_malloc( 2 * Total_size * sizeof(float_kind));
  N12  = my_malloc( 2 * Total_size * sizeof(float_kind));
  N13  = my_malloc( 2 * Total_size * sizeof(float_kind));
  FN11 = (complex_kind*) N11;
  FN12 = (complex_kind*) N12;
  FN13 = (complex_kind*) N13;
  p11  = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, FN11, N11, MPI_COMM_WORLD, FFTW_ESTIMATE);
  p12  = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, FN12, N12, MPI_COMM_WORLD, FFTW_ESTIMATE);
  p13  = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, FN13, N13, MPI_COMM_WORLD, FFTW_ESTIMATE);
#endif

  //=================================================================
  // This returns N11,N12,N13 which hold the components of
  // the vector (grad grad^{-2} density) on a grid.
  //=================================================================
  if (ThisTask == 0) printf("Calculating forces...\n");
  Forces();

#ifdef MEMORY_MODE
  my_free(density);
  my_fftw_destroy_plan(plan);
  if(modified_gravity_active) FreeMGArrays();
  for (int j = 0; j < 3; j++) Disp[j] = my_malloc(NumPart * sizeof(float));
#else
  for (int j = 0; j < 3; j++) Disp[j] = my_malloc(NumPart * sizeof(float_kind));
#endif

  //=================================================================
  // Now find the accelerations at the particle positions using 3-linear interpolation. 
  //=================================================================
  if (ThisTask == 0) printf("Calculating accelerations...\n");
  MtoParticles();

#ifdef MEMORY_MODE
  my_free(N11);
  my_free(N12);
  my_free(N13);  
  my_fftw_destroy_plan(p11);
  my_fftw_destroy_plan(p12);
  my_fftw_destroy_plan(p13);
#endif
}

//==============================================================================================
// A routine to check whether all the particles are on the correct processor and move them if not.
//==============================================================================================
void MoveParticles(void) {
  timer_start(_MoveParticles);

  //==============================================================================================
  // Note that there are some subtleties in this routine that deal with the fact in some instances there
  // may be no particles on the last N tasks depending on how the work is partioned, hence we need to 
  // skip over these tasks and copy to the correct ones. We include subtleties that deal with the fact that
  // a task may have no particles by skipping over them from the other tasks perspective and
  // setting any sendrecv commands on these tasks to null
  //==============================================================================================

  int X;
  int j;
  int neighbour, neighbour_left, neighbour_right, neighbour_count = 0;
  int send_count_max = (int)(ceil(Local_np*Nsample*Nsample*(Buffer-1.0)));
  int send_count_left = 0, send_count_right = 0;
  int recv_count_left = 0, recv_count_right = 0;
  int procdiff_left, procdiff_right, procdiffmax = 1, procdiffmaxglob = 1;
  unsigned int i;
  double scaleBox=(double)Nmesh/Box;

  //==============================================================================================
  // We assume that at least one send is needed and calculate the true number of sends needed in the first iteration.
  // (Yes, i know we shouldn't really modify the iteration counter inside the loop but it creates a good algorithm both here and
  // when we assign the particles to be copied) 
  //==============================================================================================
  for (j = 1; j <= procdiffmaxglob; j++) {

    //==============================================================================================
    // Allocate memory to hold the particles to be transfered. We assume a maximum of Local_np*Nsample*Nsample*(buffer-1.0).
    //==============================================================================================
    struct part_data * P_send_left  = (struct part_data *)my_malloc(send_count_max*sizeof(struct part_data));
    struct part_data * P_send_right = (struct part_data *)my_malloc(send_count_max*sizeof(struct part_data));

    //==============================================================================================
    // The main purpose here is to calculate how many sendrecvs we need to perform (i.e., the maximum number 
    // of tasks a particle has moved across). However, we also assume that at least one send is needed 
    // and so set up the particles to be transferred to the neighbouring tasks
    //==============================================================================================
    send_count_left = 0; send_count_right = 0;
    recv_count_left = 0; recv_count_right = 0;
    if (j <= procdiffmax) {
      for (i = 0; i < NumPart; i++) {
        X = (int)(P[i].Pos[0]*scaleBox);
        procdiff_left = 0; procdiff_right = 0;
        if (Slab_to_task[X] != ThisTask) {
          neighbour = ThisTask;
          do {
            procdiff_left++;
            neighbour--;
            if (neighbour < 0) neighbour += NTask;
            if (Local_np_table[neighbour] == 0) procdiff_left--;
          } while(Slab_to_task[X] != neighbour);
          neighbour = ThisTask;
          do {
            procdiff_right++;
            neighbour++;
            if (neighbour >= NTask) neighbour -= NTask;
            if (Local_np_table[neighbour] == 0) procdiff_right--;
          } while(Slab_to_task[X] != neighbour);
          if ((procdiff_left != 0) || (procdiff_right != 0)) {
            if (procdiff_left <= procdiff_right) {
              if (j == 1) {
                if (procdiff_left > procdiffmax) procdiffmax = procdiff_left;
              }
              if (procdiff_left == j) {
                P_send_left[send_count_left] = P[i];
                P[i] = P[NumPart-1];
                i--; NumPart--;
                send_count_left++;
                if (send_count_left >= send_count_max) {
                  printf("\nERROR: Number of particles to be sent left [%i] on task %d is greater than send_count_max [%i]\n", 
                      send_count_left, ThisTask, send_count_max);
                  printf("       You must increase the size of the buffer region.\n\n");
                 FatalError((char *)"auxPM.c Number of particles to be sent left > send_count_max");
                }
              }
            } else {
              if (j == 1) {
                if (procdiff_right > procdiffmax) procdiffmax = procdiff_right;
              }
              if (procdiff_right == j) {
                P_send_right[send_count_right] = P[i];
                P[i] = P[NumPart-1];
                i--; NumPart--;
                send_count_right++;
                if (send_count_right >= send_count_max) {
                  printf("\nERROR: Number of particles to be sent right [%i] on task %d is greater than send_count_max [%i]\n", 
                      send_count_right, ThisTask, send_count_max);
                  printf("       You must increase the size of the buffer region.\n\n");
                  FatalError((char *)"auxPM.c Number of particles to be sent right > send_count_max");
                }
              }
            }
          }
        }
      }
    } 

    //==============================================================================================
    // If we have to send to non-adjoining tasks then we have to recompute the neighbour's task number. For adjoining tasks 
    // we have already got these in the variables LeftTask and RightTask which are also used elsewhere
    //==============================================================================================
    if (j == 1) {
      neighbour_left = LeftTask;
      neighbour_right = RightTask;      
      ierr = MPI_Allreduce(&procdiffmax, &procdiffmaxglob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (ThisTask == 0) printf("Need to transfer particles %d times...\n", procdiffmaxglob);
    } else {
      if (Local_np == 0) {
        neighbour_left = MPI_PROC_NULL;
        neighbour_right = MPI_PROC_NULL;
      } else {

        neighbour_count = 0;
        neighbour_left = ThisTask;
        do {
          neighbour_left--;
          neighbour_count++;
          if(neighbour_left < 0) neighbour_left += NTask;
          if(Local_np_table[neighbour_left] == 0) neighbour_count--;
        } while(neighbour_count != j);

        neighbour_count = 0;
        neighbour_right = ThisTask;
        do {
          neighbour_right++;
          neighbour_count++;
          if(neighbour_right >= NTask) neighbour_right -= NTask;
          if(Local_np_table[neighbour_right] == 0) neighbour_count--;
        } while(neighbour_count != j);
      }
    }

    ierr = MPI_Sendrecv(&send_count_left, 1,MPI_INT,neighbour_left, 0,
                        &recv_count_right,1,MPI_INT,neighbour_right,0,
                        MPI_COMM_WORLD,&status);

    ierr = MPI_Sendrecv(&send_count_right,1,MPI_INT,neighbour_right,0,
                        &recv_count_left, 1,MPI_INT,neighbour_left, 0,
                        MPI_COMM_WORLD,&status);

    if (NumPart+recv_count_left+recv_count_right > Local_np*Nsample*Nsample*Buffer) {
      printf("\nERROR: Number of particles to be recieved on task %d is greater than available space\n", ThisTask);
      printf("       You must increase the size of the buffer region.\n\n");
      FatalError((char *)"auxPM.c Number of particles to be recieved > available space. Increase buffer");
    }

    //==============================================================================================
    // Copy across the new particles and store them at the end (of the memory). Then modify NumPart to include them.
    //==============================================================================================
    ierr = MPI_Sendrecv(&(P_send_left[0]),  send_count_left,  PartDataMPIType, neighbour_left,  0,
                        &(P[NumPart]),      recv_count_right, PartDataMPIType, neighbour_right, 0,
                        MPI_COMM_WORLD, &status);
    NumPart += recv_count_right;

    ierr = MPI_Sendrecv(&(P_send_right[0]), send_count_right, PartDataMPIType, neighbour_right, 0,
                        &(P[NumPart]),      recv_count_left,  PartDataMPIType, neighbour_left,  0,
                        MPI_COMM_WORLD, &status);
    NumPart += recv_count_left;

    my_free(P_send_left);
    my_free(P_send_right);
  }
  
  timer_stop(_MoveParticles);
  return;  
}

//==============================
// Does Cloud-in-Cell assignment.
//==============================
void PtoMesh(void) {
  timer_start(_PtoMesh);
  unsigned int i;
  unsigned int IX, IY, IZ;
  unsigned int IXneigh, IYneigh, IZneigh;
  double X, Y, Z;
  double TX, TY, TZ;
  double DX, DY, DZ;
  double scaleBox = (double)Nmesh/Box;
  double WPAR = pow((double)Nmesh / (double)Nsample,3);

  // Initialize density to -1
  for(i = 0; i < 2 * Total_size; i++) 
    density[i] = -1.0;
  
  for(i = 0; i < NumPart; i++) {

    // Scale positions to be in [0, Nmesh]
    X = P[i].Pos[0] * scaleBox;
    Y = P[i].Pos[1] * scaleBox;
    Z = P[i].Pos[2] * scaleBox;

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
    density[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]                += TX*TY*TZ;
    density[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]           += TX*TY*DZ;
    density[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]           += TX*DY*TZ;
    density[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]      += TX*DY*DZ;
    density[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]           += DX*TY*TZ;
    density[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]      += DX*TY*DZ;
    density[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]      += DX*DY*TZ;
    density[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh] += DX*DY*DZ;
  }

  //====================================================================================
  // Copy across the extra slice from the task on the left and add it to the leftmost slice
  // of the task on the right. Skip over tasks without any slices.
  //====================================================================================
  
  float_kind * temp_density = (float_kind *)my_calloc(2*alloc_slice,sizeof(float_kind));
  ierr = MPI_Sendrecv(&(density[2*last_slice]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask,0,
      &(temp_density[0]),2*alloc_slice*sizeof(float_kind),MPI_BYTE,LeftTask,0,MPI_COMM_WORLD,&status);
  if (NumPart != 0) {
    for (i = 0; i < 2 * alloc_slice; i++) density[i] += (temp_density[i] + 1.0);
  }
  my_free(temp_density);

  //====================================================================================
  // If modified gravity is active We take a copy of the density array 
  // needed to compute the fifth-force
  //====================================================================================
  if(modified_gravity_active) CopyDensityArray();

  // FFT the density field
  my_fftw_execute(plan);

#ifdef COMPUTE_POFK
  
  // Compute matter P(k) every time-step
  if(pofk_compute_every_step){
    compute_power_spectrum(P3D, aexp_global, "CDM");
  }

  // Compute matter RSD P0(k),P2(k),P4(k) every time-step
  // If we output to file in this step then we compute it then instead as this will ensure 
  // velocities are synchronized to the positions (which is not the case here!)
  int we_will_output_this_step = (timeStep_global == 0) && (i != NoutputStart_global);
  if(pofk_compute_rsd_pofk == 1 && !we_will_output_this_step){
    compute_RSD_powerspectrum(aexp_global, 0);
  }
#endif

#ifdef MASSIVE_NEUTRINOS
  // Add massive neutrinos to the total density field
  if(nu_include_massive_neutrinos){

    // When OmegaNu != 0 we need to rescale the density by OmegaCDM+Baryons/OmegaMtot
    double cdmfac = (Omega - OmegaNu) / Omega;

    // f_nu times FFT normalization factor
    double nufac_tmp   = OmegaNu / Omega * (double)(Nmesh * Nmesh * Nmesh);

    for (unsigned int i = 0; i < Local_nx; i++) {
      int iglobal = i + Local_x_start;
      for (unsigned int j = 0; j < (unsigned int)(Nmesh/2+1); j++) {
        int kmin = ((iglobal == 0) && (j == 0)) ? 1 : 0;
        for (unsigned int k = kmin; k < (unsigned int)(Nmesh/2+1); k++) {
          unsigned int coord = (i*Nmesh+j)*(Nmesh/2+1)+k;

          double dd[3], kmag;
          dd[0] = iglobal > Nmesh/2 ? iglobal-Nmesh : iglobal;
          dd[1] = j;
          dd[2] = k;
          kmag = 2.0 * PI / Box * sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2]);

          // Transform from delta_cdm(k,z=0) to delta_nu(k,z)
          double nufac = nufac_tmp * get_nu_transfer_function(kmag, aexp_global) / get_cdm_baryon_transfer_function(kmag, 1.0);
  
          P3D[coord][0] = cdmfac * P3D[coord][0] + nufac * cdelta_cdm[coord][0];
          P3D[coord][1] = cdmfac * P3D[coord][1] + nufac * cdelta_cdm[coord][1];

          if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
            coord = (i*Nmesh+(Nmesh-j))*(Nmesh/2+1)+k;
            P3D[coord][0] = cdmfac * P3D[coord][0] + nufac * cdelta_cdm[coord][0];
            P3D[coord][1] = cdmfac * P3D[coord][1] + nufac * cdelta_cdm[coord][1];
          }
        }
      }
    }
  }

#ifdef COMPUTE_POFK
  // Compute total matter P(k) every time-step
  if(pofk_compute_every_step)
    compute_power_spectrum(P3D, aexp_global, "total");
#endif

#endif

  timer_stop(_PtoMesh);
  return;
}

//===========================================
// Calculate the force grids from the density.
//===========================================
void Forces(void) {
  timer_start(_Forces);
  int dd[3];
  double RK, KK, grid_corr;
  double Scale = 2. * M_PI / Box;
  complex_kind dens;

  //==========================================================
  // Add fifth-force potential to the newtonian potential
  // For GeffG(a) models (TimeDepGeffModels) we don't allocate
  // MG arrays to save memory and we have already
  // multiplied P3D by GeffG(a) so nothing to do here
  //==========================================================
  if(modified_gravity_active && allocate_mg_arrays){
    for(unsigned int i = 0; i < Total_size; i++){
      P3D[i][0] += P3D_mgarray_two[i][0];
      P3D[i][1] += P3D_mgarray_two[i][1];
    }
  }

  //==========================================================
  // We need global values for i as opposed to local values
  // Same goes for anything that relies on i (such as RK). 
  //==========================================================
  for (unsigned int i = 0; i < Local_nx; i++) {
    int iglobal = i + Local_x_start;
    for (unsigned int j = 0; j < (unsigned int)(Nmesh/2+1); j++) {
      int kmin = 0;
      if ((iglobal == 0) && (j == 0)) {
        // Set the k = (0,0,0) mode
        FN11[0][0] = 0.0; FN11[0][1] = 0.0;
        FN12[0][0] = 0.0; FN12[0][1] = 0.0;
        FN13[0][0] = 0.0; FN13[0][1] = 0.0;
        kmin = 1;
      }
      for (unsigned int k = kmin; k < (unsigned int)(Nmesh/2+1); k++) {

        //==========================================================
        // Compute k_vec = (kx,ky,kz) and RK = |k_vec|^2
        //==========================================================
        unsigned int ind = (i*Nmesh + j)*(Nmesh/2+1) + k;
        dd[0] = iglobal > Nmesh/2 ? iglobal-Nmesh : iglobal;
        dd[1] = j;
        dd[2] = k;
        RK    = dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2];
        KK    = -1.0/RK;

        //==========================================================
        // Deconvolve the CIC window function twice (once for density, once for force interpolation)
        // and add gaussian smoothing if requested
        //==========================================================
        //grid_corr = 1.0;
        //for(int axes = 0; axes < 3; axes++)
        //  if(dd[axes] != 0) grid_corr *= sin( (PI*dd[axes]) / (double)Nmesh )/( (PI*dd[axes]) / (double)Nmesh);
        //grid_corr = pow(1.0 / grid_corr, 4.0);
        grid_corr = 1.0;

        //==========================================================
        // Gravitational potential
        //==========================================================
        dens[0] = (     P3D[ind][0]*KK*grid_corr)/pow((double)Nmesh,3);
        dens[1] = (-1.0*P3D[ind][1]*KK*grid_corr)/pow((double)Nmesh,3);

        //==========================================================
        // dens now holds the total potential so we can solve for the force. 
        //==========================================================
        FN11[ind][0] = dens[1] * dd[0] / Scale;
        FN11[ind][1] = dens[0] * dd[0] / Scale;
        FN12[ind][0] = dens[1] * dd[1] / Scale;
        FN12[ind][1] = dens[0] * dd[1] / Scale;
        FN13[ind][0] = dens[1] * dd[2] / Scale;
        FN13[ind][1] = dens[0] * dd[2] / Scale;

        //==========================================================
        // Do the mirror force along the y axis
        //==========================================================
        if ((j != (unsigned int)(Nmesh/2)) && (j != 0)) {
          int ind = (i*Nmesh + (Nmesh-j))*(Nmesh/2+1) + k;
          dd[1] = -j;

          //==========================================================
          // Gravitiational potential
          //==========================================================
          dens[0] = (     P3D[ind][0]*KK*grid_corr)/pow((double)Nmesh,3) ;
          dens[1] = (-1.0*P3D[ind][1]*KK*grid_corr)/pow((double)Nmesh,3) ;

          //==========================================================
          // dens now holds the total potential so we can solve for the force. 
          //==========================================================
          FN11[ind][0] = dens[1] * dd[0] / Scale;
          FN11[ind][1] = dens[0] * dd[0] / Scale;
          FN12[ind][0] = dens[1] * dd[1] / Scale;
          FN12[ind][1] = dens[0] * dd[1] / Scale;
          FN13[ind][0] = dens[1] * dd[2] / Scale;
          FN13[ind][1] = dens[0] * dd[2] / Scale;
        }
      }
    }
  }

  // Perform FFTs
  my_fftw_execute(p11);
  my_fftw_execute(p12);
  my_fftw_execute(p13);

  //============================================================================
  // Copy across the extra slice from the process on the right and save it at the 
  // end of the force array. Skip over tasks without any slices.
  //============================================================================
  ierr = MPI_Sendrecv(&(N11[0]), 2*alloc_slice*sizeof(float_kind), MPI_BYTE, LeftTask, 0,
      &(N11[2*last_slice]), 2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask, 0, MPI_COMM_WORLD, &status);
  ierr = MPI_Sendrecv(&(N12[0]), 2*alloc_slice*sizeof(float_kind), MPI_BYTE, LeftTask, 0,
      &(N12[2*last_slice]), 2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask, 0, MPI_COMM_WORLD, &status);
  ierr = MPI_Sendrecv(&(N13[0]), 2*alloc_slice*sizeof(float_kind), MPI_BYTE, LeftTask, 0,
      &(N13[2*last_slice]), 2*alloc_slice*sizeof(float_kind),MPI_BYTE,RightTask, 0, MPI_COMM_WORLD, &status);

  timer_stop(_Forces);
  return;
}

//===========================
// Does 3-linear interpolation
//===========================
void MtoParticles(void) {
  timer_start(_MtoParticles);
  unsigned int i;
  unsigned int IX,IY,IZ;
  unsigned int IXneigh,IYneigh,IZneigh;
  double X,Y,Z;
  double TX,TY,TZ;
  double DX,DY,DZ;
  double scaleBox = (double)Nmesh/Box;
  double WPAR = 1;

  for(int axes = 0; axes < 3; axes++)
    sumDxyz[axes] = 0;

  for(i = 0; i < NumPart; i++) {

    X = P[i].Pos[0] * scaleBox;
    Y = P[i].Pos[1] * scaleBox;
    Z = P[i].Pos[2] * scaleBox;

    IX = (unsigned int) X;
    IY = (unsigned int) Y;
    IZ = (unsigned int) Z;
   
    DX = X - (double) IX;
    DY = Y - (double) IY;
    DZ = Z - (double) IZ;
    
    TX = 1.0 - DX;
    TY = 1.0 - DY;
    TZ = 1.0 - DZ;

    DY *= WPAR;
    TY *= WPAR;

    IX -= Local_x_start;
    if(IY >= (unsigned int)Nmesh) IY = 0;
    if(IZ >= (unsigned int)Nmesh) IZ = 0;

    IXneigh = IX + 1;
    IYneigh = IY + 1;
    IZneigh = IZ + 1;
    if(IYneigh >= (unsigned int)Nmesh) IYneigh = 0;
    if(IZneigh >= (unsigned int)Nmesh) IZneigh = 0;

    Disp[0][i] = N11[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]               *TX*TY*TZ +
                 N11[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]          *TX*TY*DZ +
                 N11[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]          *TX*DY*TZ +
                 N11[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]     *TX*DY*DZ +
                 N11[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *DX*TY*TZ +
                 N11[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *DX*TY*DZ +
                 N11[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *DX*DY*TZ +
                 N11[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*DX*DY*DZ;

    Disp[1][i] = N12[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]               *TX*TY*TZ +
                 N12[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]          *TX*TY*DZ +
                 N12[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]          *TX*DY*TZ +
                 N12[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]     *TX*DY*DZ +
                 N12[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *DX*TY*TZ +
                 N12[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *DX*TY*DZ +
                 N12[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *DX*DY*TZ +
                 N12[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*DX*DY*DZ;

    Disp[2][i] = N13[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZ]               *TX*TY*TZ +
                 N13[(IX*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]          *TX*TY*DZ +
                 N13[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]          *TX*DY*TZ +
                 N13[(IX*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]     *TX*DY*DZ +
                 N13[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZ]          *DX*TY*TZ +
                 N13[(IXneigh*Nmesh+IY)*2*(Nmesh/2+1)+IZneigh]     *DX*TY*DZ +
                 N13[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZ]     *DX*DY*TZ +
                 N13[(IXneigh*Nmesh+IYneigh)*2*(Nmesh/2+1)+IZneigh]*DX*DY*DZ;

    for(int axes = 0; axes < 3; axes++)
      sumDxyz[axes] += Disp[axes][i];
  }

  // Make sumDx, sumDy and sumDz global averages
  for(int axes = 0; axes < 3; axes++){
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumDxyz[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumDxyz[axes] /= (double) TotNumPart;
  }

  timer_stop(_MtoParticles);
  return;      
}

//===============================
// Wrap the particles periodically
//===============================
#if (MEMORY_MODE || SINGLE_PRECISION)
float periodic_wrap(float x){
  while(x >= (float)Box) x -= (float)Box;
  while(x < 0) x += (float)Box;
  if (x == (float)Box) x = 0.0;
  return x;
}
#else
double periodic_wrap(double x){
  while(x >= Box) x -= Box;
  while(x < 0) x += Box;
  if (x == Box) x = 0.0;
  return x;
}
#endif

//===============
// Error message
//===============
void FatalError(char* errmsg) {
  printf("Fatal Error: [%s]\n", errmsg);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

//===========================================================================
// This catches I/O errors occuring for fwrite(). In this case we better stop.
//===========================================================================
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream) {
  size_t nwritten;
  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb) {
    printf("\nERROR: I/O error (fwrite) on task=%d has occured.\n\n", ThisTask);
    fflush(stdout);
    FatalError((char *)"auxPM.c my_fwrite | can't open file");
  }
  return nwritten;
}

//===============================================
// Create a MPI Type that described how the
// part_data struct is laid out in memory
// so we can communicate it
//===============================================

void create_MPI_type_for_Particles(MPI_Datatype *mpi_particle_type){
  const int nmax = 100;
  MPI_Aint offsets[nmax];
  MPI_Datatype types[nmax];
  int blocklengths[nmax];

  MPI_Datatype PARTICLE_FLOAT_KIND;
#ifdef MEMORY_MODE
  MPI_Type_match_size( MPI_TYPECLASS_REAL, sizeof(float), &PARTICLE_FLOAT_KIND);
#else
  MPI_Type_match_size( MPI_TYPECLASS_REAL, sizeof(float_kind), &PARTICLE_FLOAT_KIND );
#endif

  int nt = 0, nb = 0, no = 0;
#ifdef PARTICLE_ID
  offsets[no++] = offsetof(struct part_data,ID);
  types[nt++] = MPI_UNSIGNED_LONG_LONG;
  blocklengths[nb++] = 1;
#endif
  offsets[no++] = offsetof(struct part_data,Pos);
  types[nt++] = PARTICLE_FLOAT_KIND;
  blocklengths[nb++] = 3;

  offsets[no++] = offsetof(struct part_data,Vel);
  types[nt++] = PARTICLE_FLOAT_KIND;
  blocklengths[nb++] = 3;
  
  offsets[no++] = offsetof(struct part_data,D);
  types[nt++] = PARTICLE_FLOAT_KIND;
  blocklengths[nb++] = 3;
  
  offsets[no++] = offsetof(struct part_data,D2);
  types[nt++] = PARTICLE_FLOAT_KIND;
  blocklengths[nb++] = 3;

#ifdef SCALEDEPENDENT
  offsets[no++] = offsetof(struct part_data,dDdy);
  types[nt++] = PARTICLE_FLOAT_KIND;
  blocklengths[nb++] = 3;
  
  offsets[no++] = offsetof(struct part_data,dD2dy);
  types[nt++] = PARTICLE_FLOAT_KIND; 
  blocklengths[nb++] = 3;
  
  offsets[no++] = offsetof(struct part_data,coord_q);
  types[nt++] = MPI_UNSIGNED;
  blocklengths[nb++] = 1;
  
  offsets[no++] = offsetof(struct part_data,init_cpu_id);
  types[nt++] = MPI_UNSIGNED;
  blocklengths[nb++] = 1;
#endif

  if(nb != no || nb != nt || nb > nmax){
    printf("Error in create MPI struct! Numbers don't add up or increase nmax\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  MPI_Type_create_struct(nb, blocklengths, offsets, types, mpi_particle_type);
  MPI_Type_commit(mpi_particle_type);
}

void write_grid_to_file(float_kind *grid, double a, char *prefix){
  double Z;
  int Zint, Zfrac, dummy;
  char buf[300];

  // Make filename
  Z     = 1.0/a - 1.0;
  Zint  = (int)floor(Z);
  Zfrac = (int)((Z - Zint)*1000);

  // Open file
  sprintf(buf, "%s/%s_%s_z%dp%03d.%d", OutputDir, prefix, FileBase, Zint, Zfrac, ThisTask);
  FILE *fp = fopen(buf, "w");

  // Write Nmesh
  int n = (int) Nmesh;
  dummy = sizeof(int);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);
  my_fwrite(&n, sizeof(int), 1, fp);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);

  // Write Local_nx
  int nx = (int) Local_nx;
  dummy = sizeof(int);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);
  my_fwrite(&nx, sizeof(int), 1, fp);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);

  // Write Local_xstart
  int xstart = (int) Local_x_start;
  dummy = sizeof(int);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);
  my_fwrite(&xstart, sizeof(int), 1, fp);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);

  // Write the type-size
  int fsize = sizeof(float_kind);
  dummy = sizeof(int);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);
  my_fwrite(&fsize, sizeof(int), 1, fp);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);

  // Write number of grid-elements
  int npts = 2 * Total_size;
  dummy = sizeof(int);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);
  my_fwrite(&npts, sizeof(int), 1, fp);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);

  // Write grid
  dummy = npts * sizeof(float_kind);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);
  my_fwrite(grid, sizeof(float_kind), npts, fp);
  my_fwrite(&dummy, sizeof(dummy), 1, fp);

  fclose(fp);
}

