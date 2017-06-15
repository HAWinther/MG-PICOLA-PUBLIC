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

//=======================================================================================//
// This file contains most of the routines for calculating the ZA and 2LPT displacements //
//=======================================================================================//

#include "vars.h"
#include "proto.h"
#include <limits.h>

//================
// Set some units
//================
void set_units(void) {

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  G             = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble        = HUBBLE * UnitTime_in_s;

  return;
}

//==============================================
// Set up the sizes for arrays and fft routines
//==============================================
void initialize_ffts(void) {
  int *Slab_to_task_local;

  alloc_local = my_fftw_mpi_local_size_3d(Nmesh, Nmesh, Nmesh/2+1, MPI_COMM_WORLD, &Local_nx, &Local_x_start);

  Local_nx_table = my_malloc(sizeof(int) * NTask);

  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);

  if(ThisTask == 0) {
    printf("\nLocal nx\n---------------------\n");
    for(int i = 0; i < NTask; i++) printf("Task = %d: Local_nx = %d\n", i, Local_nx_table[i]);
    printf("---------------------\n");
    fflush(stdout);
  }

  //==============================================
  // Set the neighbouring tasks
  //==============================================
  if (Local_nx == 0) {
    LeftTask = MPI_PROC_NULL;
    RightTask = MPI_PROC_NULL;
  } else {
    LeftTask = ThisTask;
    do {
      LeftTask--;
      if(LeftTask < 0) LeftTask = NTask - 1;
    } while(Local_nx_table[LeftTask] == 0);

    RightTask = ThisTask;
    do {
      RightTask++;
      if(RightTask >= NTask) RightTask = 0;
    } while(Local_nx_table[RightTask] == 0);
  }

  //==============================================
  // Let each processor know which parts of the fourier grids they all have
  //==============================================
  Slab_to_task       = my_malloc(sizeof(int) * Nmesh);
  Slab_to_task_local = my_malloc(sizeof(int) * Nmesh);

  for(int i = 0; i < Nmesh; i++)    Slab_to_task_local[i] = 0;
  for(int i = 0; i < Local_nx; i++) Slab_to_task_local[Local_x_start + i] = ThisTask;

  MPI_Allreduce(Slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  //==============================================
  // Add an additional plane
  //==============================================
  alloc_slice = Nmesh*(Nmesh/2+1);
  last_slice = Local_nx*alloc_slice;
  Total_size = alloc_local+alloc_slice;

  // Sanity check
  unsigned long long max_unsigned_int = (unsigned long long) UINT_MAX;
  unsigned long long max_coord = (unsigned long long)(2 * Total_size);
  if(max_coord > max_unsigned_int){
    printf("Error: the number of gridnodes on a single task is larger than unsigned int can hold %llu > %llu !\n", max_coord, max_unsigned_int);
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  my_free(Slab_to_task_local);

  return;
}

//==============================================
// Work out which particles each Task will have
//==============================================
void initialize_parts(void) {

  int slab;
  int * Part_to_task_local;

  Local_np = 0;
  Local_p_start = Nsample;
  for (int i = 0; i < Nsample; i++) {
    slab = (int)((double)(i*Nmesh)/(double)Nsample);
    if (Slab_to_task[slab] == ThisTask) {
      Local_np++;
      if (i < Local_p_start) Local_p_start = i;
    }
  }
  
  // Sanity check
  unsigned long long max_unsigned_int = (unsigned long long) UINT_MAX;
  unsigned long long max_part_coord = (unsigned long long) Local_np * (unsigned long long) Nsample * (unsigned long long) Nsample;
  if(max_part_coord > max_unsigned_int){
    printf("Error: the number of particles on a single task is larger than unsigned int can hold %llu > %llu !\n", max_part_coord, max_unsigned_int);
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  Local_np_table = my_malloc(sizeof(int) * NTask);

  MPI_Allgather(&Local_np, 1, MPI_INT, Local_np_table, 1, MPI_INT, MPI_COMM_WORLD);

  NumPart = Local_np*Nsample*Nsample;
  TotNumPart = ((unsigned long long) Nsample) * ((unsigned long long) Nsample) *  ((unsigned long long) Nsample);

  if(ThisTask == 0) {
    printf("\nParticles\n---------------------\n");
    for(int i = 0; i < NTask; i++) printf("Task = %d: Particles = %u\n", i, (unsigned int)Local_np_table[i]*(unsigned int)Nsample*(unsigned int)Nsample);
    printf("\n----------------------\n");
    printf("Total number of particles = %llu\n\n", TotNumPart);
    fflush(stdout);
  }

  NTaskWithN = 0;
  for(int i = 0; i < NTask; i++) {
    if(Local_np_table[i] > 0) NTaskWithN++;
  }

  //==========================================================================
  // Let each processor know which parts of the particle grids they all have
  //==========================================================================
  Part_to_task       = my_malloc(sizeof(int) * Nsample);
  Part_to_task_local = my_malloc(sizeof(int) * Nsample);

  for(int i = 0; i < Nsample; i++)  Part_to_task_local[i] = 0;
  for(int i = 0; i < Local_np; i++) Part_to_task_local[Local_p_start + i] = ThisTask;

  MPI_Allreduce(Part_to_task_local, Part_to_task, Nsample, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  my_free(Part_to_task_local);

  return;
}

//======================================================================================================================
// This is the largest routine in the code and is used to generate the 2LPT initial conditions. 
// It is adapted from the 2LPTic code provided by Marc. 
//
// NOTE: We ALWAYS compute the power spectrum and displacements at redshift 0.0, regardless of 
// whether we use COLA or PM methods. These are then modified accordingly when the particle data is initialised in main.c
//======================================================================================================================
void displacement_fields(void) {
  timer_start(_DisplacementFields);
  
  if(ThisTask == 0){
    printf("\n=================================\n");
    printf("Generate displacment fields...   \n");
    printf("=================================\n");
  }

  //====================================
  // There are a LOT of parameters here:
  //====================================
  
  gsl_rng * random_generator;
  int i, j, k, n, m, p, ii, jj, kk, axes;
  unsigned int * seedtable, q, coord, bytes;
  unsigned long long nmesh3; 
  float_kind *(digrad[6]);
  float_kind *(disp[3]), *(disp2[3]);
  double u, v, w;
  double phase, ampl;
  double kvec[3], kmag, kmag2;
  double sumdis[3], sumdis2[3];
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis[3], dis2[3], maxdisp, max_disp_glob;
  complex_kind *(cdigrad[6]);
  complex_kind *(cdisp[3]), *(cdisp2[3]);
  plan_kind Forward_plan, Inverse_plan;

  //========================================================
  // General parameters for any gaussianity/non-gaussianity
  //========================================================
#ifdef GAUSSIAN
  double delta;
  double p_of_k;
#else
  float_kind *(pot);
  double t_of_k;
  double twb, phig, Beta;
  complex_kind *(cpot);
#endif

  //========================================================
  // Parameters for generic non-gaussianity
  //========================================================
#ifdef GENERIC_FNL
  int ikernel;
  double kmag0, kmagA, kmagB;
  double kerCoef, ker0, kerA, kerB;
  float_kind *(ppA), *(ppB), *(ppC);                                        
  complex_kind *(cppA), *(cppB), *(cppC);                                  
#endif

  //=============================================================
  // Parameters for equilateral or orthogonal fnl non-gaussianity
  //=============================================================
#if (EQUIL_FNL || ORTHO_FNL)
  float_kind *(partpot);
  float_kind *(p1p2p3sym), *(p1p2p3sca), *(p1p2p3nab), *(p1p2p3tre);                                                 
  complex_kind *(cpartpot);
  complex_kind *(cp1p2p3sym), *(cp1p2p3sca), *(cp1p2p3nab), *(cp1p2p3tre);
#endif                                               

  if(ThisTask == 0) {
    printf("Computing displacement fields...\n");
    fflush(stdout);
  }

  maxdisp = 0;
  
  //=============================================================
  // Initialize random number generator
  //=============================================================

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, Seed);
  if(!(seedtable = my_malloc(Nmesh * Nmesh * sizeof(unsigned int)))) FatalError((char *)"Seedtable 2LPT.c");
  for(i = 0; i < Nmesh / 2; i++) {
    for(j = 0; j < i; j++)     seedtable[i * Nmesh + j] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[j * Nmesh + i] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i; j++)     seedtable[(Nmesh - 1 - i) * Nmesh + j] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[(Nmesh - 1 - j) * Nmesh + i] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i; j++)     seedtable[i * Nmesh + (Nmesh - 1 - j)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[j * Nmesh + (Nmesh - 1 - i)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i; j++)     seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
    for(j = 0; j < i + 1; j++) seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = (unsigned int)(0x7fffffff * gsl_rng_uniform(random_generator));
  }

  //=============================================================
  // Gaussian initial conditions
  //=============================================================
#ifdef GAUSSIAN

  if(ThisTask == 0) {
    printf("Starting gaussian calculations...\n");
    fflush(stdout);
  }

  // Allocate memory and initialize arrays
  for(axes = 0, bytes = 0; axes < 3; axes++) {
    cdisp[axes] = my_malloc(sizeof(complex_kind) * Total_size);
    bytes += sizeof(complex_kind) * Total_size;
    disp[axes] = (float_kind *) cdisp[axes];
  }
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        for(axes = 0; axes < 3; axes++)  {
          coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
          cdisp[axes][coord][0] = 0.0;
          cdisp[axes][coord][1] = 0.0;
        }
      }
    }
  }

#ifdef SCALEDEPENDENT
  // Allocate grids for storing the scalar displacement fields in k-space
  cdelta_cdm  = my_malloc(sizeof(complex_kind)*Total_size);
  cdelta_cdm2 = my_malloc(sizeof(complex_kind)*Total_size);
  for(i = 0; i < Total_size; i++) {
    cdelta_cdm[i][0]  = cdelta_cdm[i][1]  = 0.0;
    cdelta_cdm2[i][0] = cdelta_cdm2[i][1] = 0.0;
  }
#endif

  if(ReadParticlesFromFile){

    //=========================================
    // Assign displacementfield from delta(k) 
    // computed from the particles we have read from file
    //=========================================

#ifdef READICFROMFILE
    AssignDisplacementField(cdisp);
#endif

  } else {

    //=====================================================================
    // If the power-spectrum we read is assumed to be for LCDM
    // then we need to rescale it. If LCDM then the ratio below is just 1
    //=====================================================================
    double sigma8_mg_over_sigma8_lcdm_pow2 = pow( mg_sigma8_enhancement(1.0), 2);

    //=========================================
    // Create IC from skratch
    //=========================================

    for(i = 0; i < Nmesh; i++) {
      ii = Nmesh - i;
      if(ii == Nmesh) ii = 0;
      if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
          (ii >= Local_x_start && ii < (Local_x_start + Local_nx))) {

        for(j = 0; j < Nmesh; j++) {
          gsl_rng_set(random_generator, seedtable[i * Nmesh + j]); 

          for(k = 0; k < Nmesh / 2; k++) {
            phase = gsl_rng_uniform(random_generator) * 2 * PI;
            do {
              ampl = gsl_rng_uniform(random_generator);
            } while(ampl == 0);

            if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2) continue;
            if(i == 0 && j == 0 && k == 0) continue;

            if(i < Nmesh / 2) {
              kvec[0] = i * 2 * PI / Box;
            } else {
              kvec[0] = -(Nmesh - i) * 2 * PI / Box;
            }

            if(j < Nmesh / 2) {
              kvec[1] = j * 2 * PI / Box;
            } else {
              kvec[1] = -(Nmesh - j) * 2 * PI / Box;
            }

            if(k < Nmesh / 2) {
              kvec[2] = k * 2 * PI / Box;
            } else {
              kvec[2] = -(Nmesh - k) * 2 * PI / Box;
            }

            kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
            kmag  = sqrt(kmag2);

            if(SphereMode == 1) {
              if(kmag * Box / (2 * PI) > Nsample / 2) continue; // select a sphere in k-space
            } else {
              if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
              if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
              if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
            }

            p_of_k  = PowerSpec(kmag);
            p_of_k *= -log(ampl);

            //====================================================================================
            // The power-spectrum we read in is assumed to be for LCDM so rescale it to get MG P(k)
            // We do this here as growth-factors are not availiable when we read in PowerSpec
            //====================================================================================
            if(modified_gravity_active && input_pofk_is_for_lcdm){
              p_of_k *= mg_pofk_ratio(kmag, 1.0);

              //=============================================
              // Since we generate at a = 1 here we need to 
              // multiply by the ratio (sigma8/sigam8_LCDM)^2
              // if the Sigma8 provided is assumed to be for
              // a corresponding LCDM simulation
              //=============================================
              if( ! input_sigma8_is_for_lcdm)
                p_of_k /= sigma8_mg_over_sigma8_lcdm_pow2;

            }

            delta = pow(Box,-1.5) * sqrt(p_of_k);  // keep at redshift 0.0
            double delta_cdm_re = delta * cos(phase);
            double delta_cdm_im = delta * sin(phase);

            if(k > 0) {
              if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                for(axes = 0; axes < 3; axes++) {
                  cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta_cdm_im;
                  cdisp[axes][coord][1] =  kvec[axes] / kmag2 * delta_cdm_re;
                }
#ifdef SCALEDEPENDENT
                cdelta_cdm[coord][0] = delta_cdm_re;
                cdelta_cdm[coord][1] = delta_cdm_im;
#endif
              }
            } else {       // k=0 plane needs special treatment
              if(i == 0) {
                if(j >= Nmesh / 2) {
                  continue;
                } else {
                  if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                    jj = Nmesh - j;  // note: j!=0 surely holds at this point
                    coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                    for(axes = 0; axes < 3; axes++) {
                      cdisp[axes][coord][0]  = -kvec[axes] / kmag2 * delta_cdm_im;
                      cdisp[axes][coord][1]  =  kvec[axes] / kmag2 * delta_cdm_re;
                    }
#ifdef SCALEDEPENDENT
                    cdelta_cdm[coord][0] = delta_cdm_re;
                    cdelta_cdm[coord][1] = delta_cdm_im;
#endif
                    // The complex conjugate here: Psi(i,j,k) = Psi*(-i,-j,-k)
                    coord = ((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
                    for(axes = 0; axes < 3; axes++) {
                      cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta_cdm_im;
                      cdisp[axes][coord][1] = -kvec[axes] / kmag2 * delta_cdm_re;
                    }
#ifdef SCALEDEPENDENT
                    cdelta_cdm[coord][0] =  delta_cdm_re;
                    cdelta_cdm[coord][1] = -delta_cdm_im;
#endif
                  }
                }
              } else {    // here comes i!=0 : conjugate can be on other processor!
                if(i >= Nmesh / 2) {
                  continue;
                } else {
                  ii = Nmesh - i;
                  jj = Nmesh - j;
                  if(ii == Nmesh) ii = 0;
                  if(jj == Nmesh) jj = 0;
                  if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                    coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                    for(axes = 0; axes < 3; axes++) {
                      cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta_cdm_im;
                      cdisp[axes][coord][1] =  kvec[axes] / kmag2 * delta_cdm_re;
                    }
#ifdef SCALEDEPENDENT
                    cdelta_cdm[coord][0] = delta_cdm_re;
                    cdelta_cdm[coord][1] = delta_cdm_im;
#endif
                  }		  
                  if(ii >= Local_x_start && ii < (Local_x_start + Local_nx)) {
                    // The complex conjugate here: Psi(i,j,k) = Psi*(-i,-j,-k)
                    coord = ((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
                    for(axes = 0; axes < 3; axes++) {
                      cdisp[axes][coord][0] = -kvec[axes] / kmag2 * delta_cdm_im;
                      cdisp[axes][coord][1] = -kvec[axes] / kmag2 * delta_cdm_re;
                    }
#ifdef SCALEDEPENDENT
                    cdelta_cdm[coord][0] =  delta_cdm_re;
                    cdelta_cdm[coord][1] = -delta_cdm_im;
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //=========================================
  // Non-gaussian initial conditions
  //=========================================
#else

  timer_start(_NonGaussianIC);
  if(ThisTask == 0) {
    printf("Starting non-gaussian calculations...\n");
    fflush(stdout);
  }

  //=========================================
  // Allocate memory and initialize arrays
  //=========================================
  bytes = 0;
  cpot = my_malloc(sizeof(complex_kind) * Total_size);
  bytes += sizeof(complex_kind) * Total_size;
  pot = (float_kind *) cpot;
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
        cpot[coord][0] = 0;
        cpot[coord][1] = 0;
      }
    }
  }

  //==========================================================
  // Ho in units of h/Mpc and c=1, i.e., internal units so far
  // Beta = 3/2 H(z)^2 a^2 Om(a) = 3/2 Ho^2 Om0 / a 
  //==========================================================
  Beta = 1.5 * Omega / FnlTime / ( INVERSE_H0_MPCH * INVERSE_H0_MPCH );

  for(i = 0; i < Nmesh; i++) {
    ii = Nmesh - i;
    if(ii == Nmesh) ii = 0;
    if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
        (ii >= Local_x_start && ii < (Local_x_start + Local_nx))) {

      for(j = 0; j < Nmesh; j++) {
        gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

        for(k = 0; k < Nmesh / 2; k++) {
          phase = gsl_rng_uniform(random_generator) * 2* PI;
          do {
            ampl = gsl_rng_uniform(random_generator);
          } while(ampl == 0);

          if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2) continue; 
          if(i == 0 && j == 0 && k == 0) continue;

          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }

          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }

          if(k < Nmesh / 2) {
            kvec[2] = k * 2 * PI / Box;
          } else {
            kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }

          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag = sqrt(kmag2);

          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue;  // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }

          phig = -log(ampl) * Anorm * pow(kmag, PrimordialIndex);            // initial normalized power
          phig = sqrt(phig) * pow(Box,-1.5) * Beta * DstartFnl / kmag2;      // amplitude of the initial gaussian potential

          if(k > 0) {
            if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
              coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
              cpot[coord][0] =  phig * sin(phase);
              cpot[coord][1] = -phig * cos(phase);
            }
          } else {  // k=0 plane needs special treatment
            if(i == 0) {
              if(j >= Nmesh / 2) {
                continue;
              } else {
                if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                  jj = Nmesh - j;   // note: j!=0 surely holds at this point

                  coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                  cpot[coord][0] =  phig * sin(phase);
                  cpot[coord][1] = - phig * cos(phase);

                  coord = ((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k; 
                  cpot[coord][0] = phig * sin(phase);
                  cpot[coord][1] = phig * cos(phase);
                }
              }
            } else {  // here comes i!=0 : conjugate can be on other processor!
              if(i >= Nmesh / 2) {
                continue;
              } else {
                ii = Nmesh - i;
                jj = Nmesh - j;
                if(ii == Nmesh) ii = 0;
                if(jj == Nmesh) jj = 0;
                if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
                  coord = ((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k;
                  cpot[coord][0] =  phig * sin(phase);
                  cpot[coord][1] = -phig * cos(phase);
                }
                if(ii >= Local_x_start && ii < (Local_x_start + Local_nx)) {
                  coord = ((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k;
                  cpot[coord][0] = phig * sin(phase);
                  cpot[coord][1] = phig * cos(phase);
                }
              }
            }
          }
        }
      }
    }
  }

  //=======================================================================================
  // For non-local models it is important to keep all factors of SQRT(-1) as done below
  // Notice also that there is a minus to convert from Bardeen to gravitational potential
  //=======================================================================================

  //==========================================================
  // Generic non-gaussian potential
  //==========================================================
#ifdef GENERIC_FNL

  // Allocate memory
  cppA = (complex_kind *) my_calloc(Total_size, sizeof(complex_kind));
  cppB = (complex_kind *) my_calloc(Total_size, sizeof(complex_kind));
  cppC = (complex_kind *) my_calloc(Total_size, sizeof(complex_kind));
  ppA  = (float_kind *) cppA;
  ppB  = (float_kind *) cppB;
  ppC  = (float_kind *) cppC;

  read_kernel_table();

  for(ikernel = 0; ikernel < NKernelTable; ikernel++) {
    kerCoef=KernelTable[ikernel].Coef;
    ker0=KernelTable[ikernel].ker0;
    kerA=KernelTable[ikernel].kerA;
    kerB=KernelTable[ikernel].kerB;

    if (ThisTask == 0) printf("\nKernel table line %d\n-------------------\n", ikernel);

    if(fabs(ker0+kerA+kerB) < 0.000000001) {
      if(ker0+kerA+kerB != 0.0) ker0 = - ( kerA + kerB );
      if(ThisTask == 0) printf("Adjusting ker0 = - (kerA + kerB), ker0=%f, kerA=%f, kerB=%f\n", ker0,kerA,kerB);
    } else {
      if(ThisTask == 0) printf("\nERROR: ker0 + kerA + kerB does not equal 0\n"); 
      FatalError((char *)"Kerneltable 2LPT.c");
    }

    if(ThisTask == 0) printf("Values: %lf %lf %lf %lf\n",kerCoef,ker0,kerA,kerB);

    for(ii = 0; ii < Local_nx; ii++) {
      for(j = 0; j < Nmesh; j++)       {
        for(k = 0; k <= Nmesh / 2 ; k++) {
          i = ii + Local_x_start;
          coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

          cppA[coord][0] = 0.0;
          cppA[coord][1] = 0.0;
          cppB[coord][0] = 0.0;
          cppB[coord][1] = 0.0;

          if((i == 0) && (j == 0) && (k == 0)) continue;

          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }

          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }

          if(k < Nmesh / 2) {
            kvec[2] = k * 2 * PI / Box;
          } else {
            kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }

          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag = sqrt(kmag2);

          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue;  // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }

          //=======================================================================
          // Eq A1 of Scoccimarro et al 1108.5512, only k dependence is relevant since 
          // normalization terms in the power cancel because ker0+kerA+kerB == 0
          //=======================================================================
          kmagA = exp((PrimordialIndex - 4.0) * kerA * log(kmag2) * 0.5); 
          kmagB = exp((PrimordialIndex - 4.0) * kerB * log(kmag2) * 0.5);

          cppA[coord][0] = kmagA * cpot[coord][0];
          cppA[coord][1] = kmagA * cpot[coord][1];

          cppB[coord][0] = kmagB * cpot[coord][0];
          cppB[coord][1] = kmagB * cpot[coord][1]; 
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(ThisTask == 0) printf("Fourier transforming initial potential ppA to configuration...\n");
    Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cppA,ppA,MPI_COMM_WORLD,FFTW_ESTIMATE);
    my_fftw_execute(Inverse_plan);
    my_fftw_destroy_plan(Inverse_plan);

    if(ThisTask == 0) printf("Fourier transforming initial potential ppB to configuration...\n");
    Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cppB,ppB,MPI_COMM_WORLD,FFTW_ESTIMATE);
    my_fftw_execute(Inverse_plan);
    my_fftw_destroy_plan(Inverse_plan);

    MPI_Barrier(MPI_COMM_WORLD);

    //===============================
    // add the terms in real space
    //===============================
    for(i = 0; i < Local_nx; i++) {
      for(j = 0; j < Nmesh; j++)    {
        for(k = 0; k < Nmesh; k++)    {
          coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
          ppA[coord] = ppA[coord]*ppB[coord];                     // ppA simply accoumulates A*B
        }
      }
    } 

    MPI_Barrier(MPI_COMM_WORLD);

    if(ThisTask == 0) printf("Fourier transforming convolution to fourier space...\n");
    Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,ppA,cppA,MPI_COMM_WORLD,FFTW_ESTIMATE);
    my_fftw_execute(Forward_plan);
    my_fftw_destroy_plan(Forward_plan);

    //=======================================================================
    // apply ker0 to the convolution of A and B
    // remove the N^3 I got by forward fourier transforming and put zero to zero mode
    //=======================================================================
    nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);    
    for(ii = 0; ii < Local_nx; ii++) {
      for(j = 0; j < Nmesh; j++)       {
        for(k = 0; k <= Nmesh / 2 ; k++) {
          i = ii + Local_x_start;
          coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

          //=============================================
          // coord = 0 is the fundamental mode in k-space
          //=============================================
          if(i == 0 && j == 0 && k == 0) {
            cppC[coord][0]=0.; 
            cppC[coord][1]=0.; 
            continue; 
          }   

          if(i < Nmesh / 2) {
            kvec[0] = i * 2 * PI / Box;
          } else {
            kvec[0] = -(Nmesh - i) * 2 * PI / Box;
          }

          if(j < Nmesh / 2) {
            kvec[1] = j * 2 * PI / Box;
          } else {
            kvec[1] = -(Nmesh - j) * 2 * PI / Box;
          }

          if(k < Nmesh / 2) {
            kvec[2] = k * 2 * PI / Box;
          } else {
            kvec[2] = -(Nmesh - k) * 2 * PI / Box;
          }

          kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          kmag = sqrt(kmag2);

          if(SphereMode == 1) {
            if(kmag * Box / (2 * PI) > Nsample / 2) continue;  // select a sphere in k-space
          } else {
            if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2) continue;
            if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2) continue;
          }

          //================================================================================
          // Eq A1 of Scoccimarro et al 1108.5512, only the k dependence is relevant since 
          // normalization terms in the power cancel because ker0+kerA+kerB == 0
          //================================================================================

          //================================================================================
          // All ikernel terms are accomulated to cppC
          //================================================================================
          kmag0 = exp((PrimordialIndex - 4.0) * ker0 * log(kmag2) * 0.5); 

          cppC[coord][0] = cppC[coord][0] + kerCoef * kmag0 * cppA[coord][0] / (double) nmesh3;   //cppA accumulated A*B earlier
          cppC[coord][1] = cppC[coord][1] + kerCoef * kmag0 * cppA[coord][1] / (double) nmesh3;   //cppA accumulated A*B earlier               
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  //==========================================
  // add the linear (i.e, non quadratic) part
  //==========================================
  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2; k++)  {
        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        //=================================================
        // cpot has linear part + cppC has quadratic part 
        //=================================================
        cpot[coord][0] = cpot[coord][0] + Fnl * cppC[coord][0];   
        cpot[coord][1] = cpot[coord][1] + Fnl * cppC[coord][1];             
      }
    }
  } 
  MPI_Barrier(MPI_COMM_WORLD);

  // Free up memory
  my_free(cppA);
  my_free(cppB);
  my_free(cppC);

  //===================================================================================
  // Specific type of non-gaussianity (local, equilateral or orthogonal) for n_s=1 only
  //===================================================================================
#else

  //===========================
  // Local primordial potential
  //===========================
#ifdef LOCAL_FNL  

  if(ThisTask == 0) printf("Fourier transforming initial potential to configuration...\n");
  Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpot,pot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Inverse_plan);
  my_fftw_destroy_plan(Inverse_plan);

  //==============================================
  // square the potential in configuration space
  //==============================================
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)    {
      for(k = 0; k < Nmesh; k++)    {
        coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
        pot[coord] = pot[coord] + Fnl * pot[coord]*pot[coord];
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0) printf("Fourier transforming squared potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,pot,cpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  //===================================================================================
  // remove the N^3 I got by forward fourier transforming and put zero to zero mode 
  //===================================================================================
  nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);    
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
        cpot[coord][0] /= (double) nmesh3; 
        cpot[coord][1] /= (double) nmesh3; 
      }
    }
  }

  if(ThisTask == 0) {
    cpot[0][1] = 0.0;
    cpot[0][0] = 0.0; 
  }

  //=================================================
  // Non-local primordial potential and initialize
  //=================================================
#else

  // Allocate memory for partpotential
  cpartpot   = my_malloc(bytes = sizeof(complex_kind) * Total_size);
  cp1p2p3sym = my_malloc(bytes = sizeof(complex_kind) * Total_size);
  cp1p2p3sca = my_malloc(bytes = sizeof(complex_kind) * Total_size);
  cp1p2p3nab = my_malloc(bytes = sizeof(complex_kind) * Total_size);
  cp1p2p3tre = my_malloc(bytes = sizeof(complex_kind) * Total_size);
  partpot    = (float_kind *) cpartpot;
  p1p2p3sym  = (float_kind *) cp1p2p3sym;
  p1p2p3sca  = (float_kind *) cp1p2p3sca;
  p1p2p3nab  = (float_kind *) cp1p2p3nab;
  p1p2p3tre  = (float_kind *) cp1p2p3tre;

  // Initialize
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
        cp1p2p3sym[coord][0] = 0;
        cp1p2p3sym[coord][1] = 0;
        cp1p2p3sca[coord][0] = 0;
        cp1p2p3sca[coord][1] = 0;
        cp1p2p3nab[coord][0] = 0;
        cp1p2p3nab[coord][1] = 0;
        cp1p2p3tre[coord][0] = 0;
        cp1p2p3tre[coord][1] = 0;
        cpartpot[coord][0] = 0;
        cpartpot[coord][1] = 0;
      }
    }
  }

  // multiply by k  
  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2; k++)  {

        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if(i < Nmesh / 2) {
          kvec[0] = i * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - i) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        kmag = sqrt(kmag2);

        cpartpot[coord][0] = kmag * cpot[coord][0];
        cpartpot[coord][1] = kmag * cpot[coord][1];

        cp1p2p3nab[coord][0] = kmag2 * cpot[coord][0];
        cp1p2p3nab[coord][1] = kmag2 * cpot[coord][1]; 
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //=================================================
  // Fourier transform back to real space
  //=================================================

  if(ThisTask == 0) printf("Fourier transforming initial potential to configuration...\n");
  Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpot,pot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Inverse_plan);
  my_fftw_destroy_plan(Inverse_plan);

  if(ThisTask == 0) printf("Fourier transforming partpotential to configuration...\n");
  Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cpartpot,partpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Inverse_plan);
  my_fftw_destroy_plan(Inverse_plan);

  if(ThisTask == 0) printf("Fourier transforming nabpotential to configuration...\n");
  Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cp1p2p3nab,p1p2p3nab,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Inverse_plan);
  my_fftw_destroy_plan(Inverse_plan);

  MPI_Barrier(MPI_COMM_WORLD);

  //=================================================
  // Multiplying terms in real space
  //=================================================
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)    {
      for(k = 0; k < Nmesh; k++)    {
        coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k; 
        p1p2p3sym[coord] = partpot[coord]*partpot[coord];
        p1p2p3sca[coord] = pot[coord]*partpot[coord];   
        p1p2p3nab[coord] = pot[coord]*p1p2p3nab[coord];
        p1p2p3tre[coord] = p1p2p3nab[coord]*partpot[coord];
        partpot[coord] = pot[coord]*pot[coord];                 // NOTE: now partpot is potential squared
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0) printf("Fourier transforming potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,pot,cpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  if(ThisTask == 0) printf("Fourier transforming squared potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,partpot,cpartpot,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  if(ThisTask == 0) printf("Fourier transforming p1p2p3sym potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3sym,cp1p2p3sym,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  if(ThisTask == 0) printf("Fourier transforming p1p2p3sca potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3sca,cp1p2p3sca,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  if(ThisTask == 0) printf("Fourier transforming p1p2p3nab potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3nab,cp1p2p3nab,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  if(ThisTask == 0) printf("Fourier transforming p1p2p3tre potential ...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,p1p2p3tre,cp1p2p3tre,MPI_COMM_WORLD,FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  MPI_Barrier(MPI_COMM_WORLD);

  //====================================================================================
  // Divide by appropiate k's, sum terms according to non-local model 
  // remove the N^3 I got by forward fourier transforming and put zero to zero mode 
  //====================================================================================
  nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);    
  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2 ; k++) {
        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if(i < Nmesh / 2) {
          kvec[0] = i * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - i) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        kmag = sqrt(kmag2);

        if(i == 0 && j == 0 && k == 0) {
          cpot[0][0]=0.;
          cpot[0][1]=0.;
          continue;
        }

        //===================     
        // F_nl equilateral
        //===================     
#ifdef EQUIL_FNL

        cpot[coord][0] = cpot[coord][0] + Fnl * (-3*cpartpot[coord][0] - 2*cp1p2p3sym[coord][0] / kmag2 + 4*cp1p2p3sca[coord][0] /kmag + 2*cp1p2p3nab[coord][0] / kmag2);
        cpot[coord][1] = cpot[coord][1] + Fnl * (-3*cpartpot[coord][1] - 2*cp1p2p3sym[coord][1] / kmag2 + 4*cp1p2p3sca[coord][1] /kmag + 2*cp1p2p3nab[coord][1] / kmag2); 
        cpot[coord][0] /= (double) nmesh3; 
        cpot[coord][1] /= (double) nmesh3; 

#endif

        //===================     
        // F_nl othogonal
        //===================     
#ifdef ORTHO_FNL 

        cpot[coord][0] = cpot[coord][0] + Fnl * (-9*cpartpot[coord][0] - 8*cp1p2p3sym[coord][0] / kmag2 + 10*cp1p2p3sca[coord][0] /kmag + 8*cp1p2p3nab[coord][0] / kmag2);
        cpot[coord][1] = cpot[coord][1] + Fnl * (-9*cpartpot[coord][1] - 8*cp1p2p3sym[coord][1] / kmag2 + 10*cp1p2p3sca[coord][1] /kmag + 8*cp1p2p3nab[coord][1] / kmag2);
        cpot[coord][0] /= (double) nmesh3; 
        cpot[coord][1] /= (double) nmesh3; 

#endif

        if(i == 0 && j == 0 && k == 0) {
          cpot[0][0]=0.;
          cpot[0][1]=0.;
          continue;
        }
      }
    }
  }

  my_free(cpartpot);
  my_free(cp1p2p3sym);
  my_free(cp1p2p3sca);
  my_free(cp1p2p3nab);
  my_free(cp1p2p3tre);

#endif

  //===================================================
  // Compute displacements using nongaussian potential
  //===================================================
#endif 

  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask==0) {
    printf("Computing displacement using non-gaussian potential...\n");
    fflush(stdout);
  }

  // Allocate memory
  for(axes = 0, bytes = 0; axes < 3; axes++) {
    cdisp[axes] = my_malloc(sizeof(complex_kind) * Total_size);
    bytes += sizeof(complex_kind) * Total_size;
    disp[axes] = (float_kind *) cdisp[axes];
  }

  if(ThisTask == 0) {
    printf("Starting axes = %d...\n", axes);
    fflush(stdout);
  }

  // First, clean the array
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        for(axes = 0; axes < 3; axes++) {
          cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][0] = 0;
          cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][1] = 0;
        }
      }
    }
  } 

  for(ii = 0; ii < Local_nx; ii++) {
    for(j = 0; j < Nmesh; j++)       {
      for(k = 0; k <= Nmesh / 2 ; k++) {
        i = ii + Local_x_start;
        coord = (ii * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if(i < Nmesh / 2) {
          kvec[0] = i * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - i) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        kmag = sqrt(kmag2);

        t_of_k = TransferFunc(kmag);
        twb = t_of_k / (Beta * DstartFnl);   

        for(axes = 0; axes < 3; axes++) {
          cdisp[axes][coord][1] = kvec[axes] * twb * cpot[coord][0];
          cdisp[axes][coord][0] = - kvec[axes] * twb * cpot[coord][1];
        }
      }
    }
  } 

  my_free(cpot);
  timer_stop(_NonGaussianIC);
#endif

  //================================================================
  // Compute the displacements (identical regardless of gaussianity)
  // At this point, cdisp[axes] contains the complex Zeldovich displacement
  //================================================================

  if(ThisTask == 0) {
    printf("Computing 2LPT displacements...\n");
    fflush(stdout);
  }

  // Allocate memory
  for(i = 0; i < 6; i++) {
    cdigrad[i] = my_malloc(sizeof(complex_kind) * Total_size);
    bytes += sizeof(complex_kind) * Total_size;
    digrad[i]  = (float_kind *) cdigrad[i];
  }

  //================================================================
  // Compute displacement gradient
  //================================================================
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if((i + Local_x_start) < Nmesh / 2) {
          kvec[0] = (i + Local_x_start) * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
        }  

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }    

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        //================================================================
        // Derivatives of ZA displacement
        // d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i
        //================================================================
        cdigrad[0][coord][0] = -cdisp[0][coord][1] * kvec[0]; // disp0,0
        cdigrad[0][coord][1] =  cdisp[0][coord][0] * kvec[0];

        cdigrad[1][coord][0] = -cdisp[0][coord][1] * kvec[1]; // disp0,1
        cdigrad[1][coord][1] =  cdisp[0][coord][0] * kvec[1];

        cdigrad[2][coord][0] = -cdisp[0][coord][1] * kvec[2]; // disp0,2
        cdigrad[2][coord][1] =  cdisp[0][coord][0] * kvec[2];

        cdigrad[3][coord][0] = -cdisp[1][coord][1] * kvec[1]; // disp1,1
        cdigrad[3][coord][1] =  cdisp[1][coord][0] * kvec[1];

        cdigrad[4][coord][0] = -cdisp[1][coord][1] * kvec[2]; // disp1,2
        cdigrad[4][coord][1] =  cdisp[1][coord][0] * kvec[2];

        cdigrad[5][coord][0] = -cdisp[2][coord][1] * kvec[2]; // disp2,2
        cdigrad[5][coord][1] =  cdisp[2][coord][0] * kvec[2];
      }
    }
  }

  if(ThisTask == 0) printf("Fourier transforming displacement gradient...\n");
  for(i = 0; i < 6; i++) {
    Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdigrad[i],digrad[i],MPI_COMM_WORLD,FFTW_MEASURE);
    my_fftw_execute(Inverse_plan);
    my_fftw_destroy_plan(Inverse_plan);
  }

  //================================================================
  // Compute second order source and store it in digrad[3]
  //================================================================
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)    {
      for(k = 0; k < Nmesh; k++)    {
        coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;
        digrad[3][coord] = digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])+digrad[3][coord]*digrad[5][coord]
          -digrad[1][coord]*digrad[1][coord]-digrad[2][coord]*digrad[2][coord]-digrad[4][coord]*digrad[4][coord];
      }
    }
  }

  if(ThisTask == 0) printf("Fourier transforming second order source...\n");
  Forward_plan = my_fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,digrad[3],cdigrad[3], MPI_COMM_WORLD, FFTW_ESTIMATE);
  my_fftw_execute(Forward_plan);
  my_fftw_destroy_plan(Forward_plan);

  //========================================================================================
  // The memory allocated for cdisp2[0], [1], and [2] will be used for 2nd order displacements
  // Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later
  //========================================================================================
  for(axes = 0; axes < 3; axes++) {
    cdisp2[axes] = cdigrad[axes]; 
    disp2[axes] = (float_kind *) cdisp2[axes];
  }
  my_free(cdigrad[4]); 
  my_free(cdigrad[5]); 

  //========================================================================================
  // Solve Poisson eq. and calculate 2nd order displacements
  //========================================================================================
  for(i = 0; i < Local_nx; i++) {
    for(j = 0; j < Nmesh; j++)     {
      for(k = 0; k <= Nmesh / 2; k++) {
        coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;

        if((i + Local_x_start) < Nmesh / 2) {
          kvec[0] = (i + Local_x_start) * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];

#ifdef SCALEDEPENDENT
        // The minus sign is to have it on the same form as the 1LPT computed above
        if(kmag2 > 0.0) {
          cdelta_cdm2[coord][0] = -cdigrad[3][coord][0];
          cdelta_cdm2[coord][1] = -cdigrad[3][coord][1];
        } else {
          cdelta_cdm2[coord][0] =  cdelta_cdm2[coord][1] = 0.0;
        }
#else
        //==================================================
        // cdisp2 = source * k / (sqrt(-1) k^2) */
        //==================================================
        for(axes = 0; axes < 3; axes++) {
          if(kmag2 > 0.0) {
            cdisp2[axes][coord][0] =  cdigrad[3][coord][1] * kvec[axes] / kmag2;
            cdisp2[axes][coord][1] = -cdigrad[3][coord][0] * kvec[axes] / kmag2;
          } else {
            cdisp2[axes][coord][0] = cdisp2[axes][coord][1] = 0.0;
          }
        }
#endif
      }
    }
  }
  
#ifdef SCALEDEPENDENT

  // We do not need to perform the last steps so just free memory and return
  my_free(cdigrad[3]);
  for(axes = 0; axes < 3; axes++) my_free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) my_free(cdisp2[axes]);
  gsl_rng_free(random_generator);
  my_free(seedtable);
  nmesh3 = ((unsigned long long ) Nmesh ) * ((unsigned long long) Nmesh) *  ((unsigned long long) Nmesh);
  timer_stop(_DisplacementFields);
  return;

#endif

  // Free cdigrad
  my_free(cdigrad[3]);

  //======================================================================
  // Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements
  //======================================================================
  for(axes = 0; axes < 3; axes++) {
    if(ThisTask == 0) printf("Fourier transforming displacements, axis %d\n",axes);
    Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp[axes],disp[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    my_fftw_execute(Inverse_plan);
    my_fftw_destroy_plan(Inverse_plan);
    Inverse_plan = my_fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp2[axes],disp2[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    my_fftw_execute(Inverse_plan);
    my_fftw_destroy_plan(Inverse_plan);

    // now get the plane on the right side from neighbour on the right and send the left plane

    // Send ZA disp
    MPI_Sendrecv(&(disp[axes][0]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,LeftTask,10,
        &(disp[axes][2*last_slice]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,RightTask,10,MPI_COMM_WORLD,&status);

    // Send 2nd order disp
    MPI_Sendrecv(&(disp2[axes][0]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,LeftTask,10,
        &(disp2[axes][2*last_slice]),sizeof(float_kind)*2*alloc_slice,MPI_BYTE,RightTask,10,MPI_COMM_WORLD,&status);     
  }

  // Read-out displacements
  gsl_rng_free(random_generator);
  my_free(seedtable);

  nmesh3 = ((unsigned long long ) Nmesh ) * ((unsigned long long) Nmesh) *  ((unsigned long long) Nmesh);

  for(axes = 0; axes < 3; axes++) {
#ifdef MEMORY_MODE
    ZA[axes]  = my_malloc(NumPart*sizeof(float));
    LPT[axes] = my_malloc(NumPart*sizeof(float));
#else
    ZA[axes]  = my_malloc(NumPart*sizeof(float_kind));
    LPT[axes] = my_malloc(NumPart*sizeof(float_kind));
#endif
    sumdis[axes] = 0;
    sumdis2[axes] = 0;    
  }

  for (n = 0; n < Local_np; n++) {
    for (m = 0; m < Nsample; m++) {
      for (p = 0; p < Nsample; p++) {
        coord = (n * Nsample + m) * (Nsample) + p;

        u = (double)((n+Local_p_start)*Nmesh)/(double)Nsample;
        v = (double)(m*Nmesh)/(double)Nsample;
        w = (double)(p*Nmesh)/(double)Nsample;

        i = (int) u;
        j = (int) v;
        k = (int) w;

        if(i == (Local_x_start + Local_nx)) i = (Local_x_start + Local_nx) - 1;
        if(i < Local_x_start)               i = Local_x_start;
        if(j == Nmesh)                      j = Nmesh - 1;
        if(k == Nmesh)                      k = Nmesh - 1;

        u -= i;
        v -= j;
        w -= k;

        i -= Local_x_start;
        ii = i + 1;
        jj = j + 1;
        kk = k + 1;

        if(jj >= Nmesh) jj -= Nmesh;
        if(kk >= Nmesh) kk -= Nmesh;

        f1 = (1 - u) * (1 - v) * (1 - w);
        f2 = (1 - u) * (1 - v) * (w);
        f3 = (1 - u) * (v) * (1 - w);
        f4 = (1 - u) * (v) * (w);
        f5 = (u) * (1 - v) * (1 - w);
        f6 = (u) * (1 - v) * (w); 
        f7 = (u) * (v) * (1 - w);
        f8 = (u) * (v) * (w);

        for(axes = 0; axes < 3; axes++) {
          dis[axes] = disp[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + k]  * f1 +
            disp[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
            disp[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + k]  * f3 +
            disp[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
            disp[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + k]  * f5 +
            disp[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
            disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k]  * f7 +
            disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

          dis2[axes] = disp2[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + k]  * f1 +
            disp2[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
            disp2[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + k]  * f3 +
            disp2[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
            disp2[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + k]  * f5 +
            disp2[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
            disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k]  * f7 +
            disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

          dis2[axes] /= (double) nmesh3; 

          ZA[axes][coord]   =  dis[axes];
          LPT[axes][coord]  =  -3./7.*dis2[axes];
          sumdis[axes]  += dis[axes];
          sumdis2[axes] += -3./7.*dis2[axes];

          if(fabs(dis[axes] - 3./7. * dis2[axes]) > maxdisp) maxdisp = fabs(dis[axes] - 3./7. * dis2[axes]);
        }
      }
    }
  }

  //======================================================================
  // Make sure the average of the displacements is zero.
  //======================================================================
  for(axes = 0; axes < 3; axes++) {
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumdis[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumdis2[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumdis[axes]  /= (double)TotNumPart;
    sumdis2[axes] /= (double)TotNumPart;
  }
  for(axes = 0; axes < 3; axes++) my_free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) my_free(cdisp2[axes]);

  for(q = 0; q < NumPart; q++) {
    for(axes = 0; axes < 3; axes++) {
      // In original code there was += instead of -= which does not make sense
      ZA[axes][q]  -= sumdis[axes];
      LPT[axes][q] -= sumdis2[axes];
    }
  }

  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (ThisTask == 0) {
    printf("Calculated Zeldovich and 2LPT displacements...\n");
    printf("Maximum displacement = %lf kpc/h (%lf in units of the particle separation)...\n\n",max_disp_glob, max_disp_glob / (Box / Nmesh));
    fflush(stdout);
  }

  timer_stop(_DisplacementFields);
  return;
}

#ifdef SCALEDEPENDENT

//==========================================================================//
//                                                                          //
//  MG-PICOLA written by Hans Winther (ICG Portsmouth) March 2017           //
//                                                                          //
// The rest of this file contains methods to compute scale-dependent        //
// displacement fields and communication between different tasks   .        //
//                                                                          //
//==========================================================================//

//====================================================================================
// This routine takes the stored displacement-field in k-space and makes it at the given redshift
// Assuming [ZA_D] has been allocated
// Assuming we already have stored the initial displacement-field in [cdisp_store]
//====================================================================================

void from_cdisp_store_to_ZA(double A, double AF, double AFF, int fieldtype, int LPTorder){
  unsigned long long nmesh3 = ((unsigned long long) Nmesh) * ((unsigned long long) Nmesh ) * ((unsigned long long) Nmesh);    

  // Growth-factor to given LPT order
  double (*func_growth_D_scaledependent)(double, double) = NULL;
  double (*func_growth_dDdy_scaledependent)(double, double) = NULL;
  double (*func_growth_ddDddy_scaledependent)(double, double) = NULL;
  double normfactor = 1.0;
  if(LPTorder == 1){
    func_growth_D_scaledependent      = &growth_D_scaledependent;
    func_growth_dDdy_scaledependent   = &growth_dDdy_scaledependent;
    func_growth_ddDddy_scaledependent = &growth_ddDddy_scaledependent;
    normfactor = 1.0; 
  } else if(LPTorder == 2){
    func_growth_D_scaledependent      = &growth_D2_scaledependent;
    func_growth_dDdy_scaledependent   = &growth_dD2dy_scaledependent;
    func_growth_ddDddy_scaledependent = &growth_ddD2ddy_scaledependent;
    normfactor = -3.0 / 7.0 / (double) nmesh3; 
  }

  if( !(LPTorder == 1 || LPTorder == 2) ){
    printf("Error in from_cdisp_store_to_ZA: LPTorder [%i] is not supported [1] or [2] are only options\n", LPTorder);
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  // Pointer to the stored [LPTorder]LPT displacment field in k-space
#ifdef SCALEDEPENDENT
  complex_kind *stored_delta;
  if(LPTorder == 1){
    stored_delta = cdelta_cdm;
  } else if(LPTorder == 2){
    stored_delta = cdelta_cdm2;
  }
#else
  complex_kind *(stored_disp_field[3]);
  if(LPTorder == 1){
    for(int axes = 0; axes < 3; axes++)
      stored_disp_field[axes] = cdisp_store[axes];
  } else if(LPTorder == 2){
    for(int axes = 0; axes < 3; axes++)
      stored_disp_field[axes] = cdisp2_store[axes];
  }
#endif

  // Allocate temporary memory
  complex_kind *(cdisp_D[3]);
  double sumdis_D[3];
  for(int axes = 0; axes < 3; axes++){
    disp_D[axes]   = my_malloc(sizeof(float_kind) * 2 * Total_size);
    cdisp_D[axes]  = (complex_kind *) disp_D[axes];
    sumdis_D[axes] = 0.0;
  }

  // Multiply by growth-factor D(k, A)
  for(int i = 0; i < Local_nx; i++) {
    for(int j = 0; j < Nmesh; j++)     {
      for(int k = 0; k <= Nmesh / 2; k++) {
        unsigned int coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;

        double kvec[3];
        if((i + Local_x_start) < Nmesh / 2) {
          kvec[0] = (i + Local_x_start) * 2 * PI / Box;
        } else {
          kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
        }

        if(j < Nmesh / 2) {
          kvec[1] = j * 2 * PI / Box;
        } else {
          kvec[1] = -(Nmesh - j) * 2 * PI / Box;
        }

        if(k < Nmesh / 2) {
          kvec[2] = k * 2 * PI / Box;
        } else {
          kvec[2] = -(Nmesh - k) * 2 * PI / Box;
        }

        // Square norm of wave-vector in units of h/Mpc
        double kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
        double kmag = sqrt(kmag2);

        // Fetch growth factors
        double growth_factor = 1.0;
        if(fieldtype == FIELD_D)      growth_factor = normfactor * func_growth_D_scaledependent(kmag, A);
        if(fieldtype == FIELD_dDdy)   growth_factor = normfactor * func_growth_dDdy_scaledependent(kmag, A);
        if(fieldtype == FIELD_ddDddy) growth_factor = normfactor * func_growth_ddDddy_scaledependent(kmag, A);
        if(fieldtype == FIELD_deltaD) growth_factor = normfactor * (func_growth_D_scaledependent(kmag, AFF) - func_growth_D_scaledependent(kmag, A));

        for(int axes = 0; axes < 3; axes++) {
          if(kmag2 > 0.0) {
#ifdef SCALEDEPENDENT
            cdisp_D[axes][coord][0] = -stored_delta[coord][1] * kvec[axes]/kmag2 * growth_factor;
            cdisp_D[axes][coord][1] =  stored_delta[coord][0] * kvec[axes]/kmag2 * growth_factor;
#else
            cdisp_D[axes][coord][0] = stored_disp_field[axes][coord][0] * growth_factor;
            cdisp_D[axes][coord][1] = stored_disp_field[axes][coord][1] * growth_factor;
#endif
          } else {
            cdisp_D[axes][coord][0] = 0.0;
            cdisp_D[axes][coord][1] = 0.0;
          }
        }
      }
    }
  }

  // Fourier transform to k-space
  for(int axes = 0; axes < 3; axes++) {
    // FFT and copy over slice
    plan_kind Inverse_plan_D = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, 
        cdisp_D[axes], disp_D[axes], MPI_COMM_WORLD, FFTW_ESTIMATE);
    my_fftw_execute(Inverse_plan_D);
    my_fftw_destroy_plan(Inverse_plan_D);
    MPI_Sendrecv(&(disp_D[axes][0]),   sizeof(float_kind) * 2 * alloc_slice, MPI_BYTE, LeftTask,  10,
        &(disp_D[axes][2*last_slice]), sizeof(float_kind) * 2 * alloc_slice, MPI_BYTE, RightTask, 10, MPI_COMM_WORLD, &status);
  }

  // Make the real-space Lagrangian displacement vectors
  char fieldname[10];
  if(fieldtype == FIELD_D)      sprintf(fieldname,"%s","[D     ]");
  if(fieldtype == FIELD_dDdy)   sprintf(fieldname,"%s","[dDdy  ]");
  if(fieldtype == FIELD_ddDddy) sprintf(fieldname,"%s","[ddDddy]");
  if(fieldtype == FIELD_deltaD) sprintf(fieldname,"%s","[deltaD]");
  if(ThisTask == 0) printf("Assigning scale-dependent %iLPT displacementfield %s to particles\n", LPTorder, fieldname);
  double maxdisp_glob = 0, maxdisp = 0;
  for (int n = 0; n < Local_np; n++) {
    for (int m = 0; m < Nsample; m++) {
      for (int p = 0; p < Nsample; p++) {
        unsigned int coord = (n * Nsample + m) * (Nsample) + p;

        // Try to understand why coord != coord when Nsample = Nmesh

        double u = (double)((n+Local_p_start)*Nmesh)/(double)Nsample;
        double v = (double)(m*Nmesh)/(double)Nsample;
        double w = (double)(p*Nmesh)/(double)Nsample;

        int i = (int) u;
        int j = (int) v;
        int k = (int) w;

        if(i == (Local_x_start + Local_nx)) i = (Local_x_start + Local_nx) - 1;
        if(i < Local_x_start)               i = Local_x_start;
        if(j == Nmesh)                      j = Nmesh - 1;
        if(k == Nmesh)                      k = Nmesh - 1;

        u -= i;
        v -= j;
        w -= k;

        i -= Local_x_start;
        int ii = i + 1;
        int jj = j + 1;
        int kk = k + 1;

        if(jj >= Nmesh) jj -= Nmesh;
        if(kk >= Nmesh) kk -= Nmesh;

        double f1 = (1 - u) * (1 - v) * (1 - w);
        double f2 = (1 - u) * (1 - v) * (w);
        double f3 = (1 - u) * (v) * (1 - w);
        double f4 = (1 - u) * (v) * (w);
        double f5 = (u) * (1 - v) * (1 - w);
        double f6 = (u) * (1 - v) * (w); 
        double f7 = (u) * (v) * (1 - w);
        double f8 = (u) * (v) * (w);

        // Trilinear interpolation
        double dis_D[3];
        for(int axes = 0; axes < 3; axes++) {

          dis_D[axes] = disp_D[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + k ] * f1 +
                        disp_D[axes][(i * Nmesh + j)   * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
                        disp_D[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + k ] * f3 +
                        disp_D[axes][(i * Nmesh + jj)  * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
                        disp_D[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + k ] * f5 +
                        disp_D[axes][(ii * Nmesh + j)  * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
                        disp_D[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k ] * f7 +
                        disp_D[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

          sumdis_D[axes]   += dis_D[axes];
          ZA_D[axes][coord] = dis_D[axes];

          if( fabs(dis_D[axes]) > maxdisp) maxdisp = dis_D[axes];
        }
      }
    }
  }

  MPI_Reduce(&maxdisp, &maxdisp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(ThisTask == 0 && fieldtype == FIELD_D)
    printf("Maximum %iLPT displacement = %lf kpc/h (%lf in units of the particle separation)\n", LPTorder, maxdisp_glob, maxdisp_glob / (Box / Nmesh));

  // Communicate sum of displacement over all particles on all CPUs
  for(int axes = 0; axes < 3; axes++) {
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumdis_D[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumdis_D[axes]  /= (double) TotNumPart;
  }

  // Make sure sum of displacements is zero
  int NumPart_init = Local_np * Nsample * Nsample;
  for(int axes = 0; axes < 3; axes++){
    for (int coord = 0; coord < NumPart_init; coord++) {
      ZA_D[axes][coord] -= sumdis_D[axes];
    }
  }

  // Free up memory
  for(int axes = 0; axes < 3; axes++){
    my_free(disp_D[axes]);
  }
}

//====================================================================
// The data we send to CPU[i] requesting displacment-fields at [coord_q]
// needed for particle [pid] (on this CPU not [i])
//====================================================================
struct disp_comm_data{
  unsigned int pid;
  unsigned int coord_q;
};

//====================================================================
// The data we get back from CPU[i]
// Use to assign displacment-fields to particle [pid] (on this CPU not [i])
//====================================================================
struct disp_comm_data_return{
  unsigned int pid;
  float_kind D[3];
};

//====================================================================
// This routine assigns the density field to the particles
//====================================================================
void assign_displacment_field_to_particles(double A, double AF, double AFF, int fieldtype, int LPTorder){
  long long allocated_bytes = 0;
  int comm_verbose = 0;

  if( !(LPTorder == 1 || LPTorder == 2) ){
    printf("Error in assign_displacment_field_to_particles: LPTorder [%i] is not supported [1] or [2] are only options\n", LPTorder);
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  // First we allocate and compute ZA_D
  for(int axes = 0; axes < 3; axes++) {
    ZA_D[axes]       = my_malloc(Nsample*Nsample*Local_np*sizeof(float_kind));
    allocated_bytes += Nsample*Nsample*Local_np*sizeof(float_kind); 
  }

  // Compute the displacment fields ZA_D, ZA_dDdy and ZA_ddDddy 
  from_cdisp_store_to_ZA(A, AF, AFF, fieldtype, LPTorder);

  // Allocate communication count buffer
  int *n_ThisTask_needs_from_cpu = my_malloc(sizeof(int) * NTask);
  for(int i = 0; i < NTask; i++) { 
    n_ThisTask_needs_from_cpu[i] = 0;
  }

  // Count how much data we need to get from other CPUs
  for(int i = 0; i < NumPart; i++){
    unsigned int cpu_id = P[i].init_cpu_id;
    if( cpu_id != ThisTask ){
      n_ThisTask_needs_from_cpu[cpu_id]++;
    }
  }

  // Allocate request_buffer[i] is list of what ThisTask requests from cpu[i]
  struct disp_comm_data **request_buffer = my_malloc(sizeof(struct disp_comm_data *) * NTask);
  for(int i = 0; i < NTask; i++) { 
    request_buffer[i] = my_malloc(sizeof(struct disp_comm_data) * n_ThisTask_needs_from_cpu[i]);
    for(int j = 0; j < n_ThisTask_needs_from_cpu[i]; j++){
      request_buffer[i][j].pid     = 0;
      request_buffer[i][j].coord_q = 0;
    }
    allocated_bytes += sizeof(struct disp_comm_data) * n_ThisTask_needs_from_cpu[i];
  }

  // Reset counter
  for(int i = 0; i < NTask; i++)
    n_ThisTask_needs_from_cpu[i] = 0;

  // Loop over all particles and make request_buffer
  // If particle is on current CPU then assign displacment-field to particle
  for(int i = 0; i < NumPart; i++){
    int cpu_id           = P[i].init_cpu_id;
    unsigned int coord_q = P[i].coord_q;

    if( cpu_id == ThisTask ){

      // Assign displacment-fields to particle
      if(LPTorder == 1){
        for(int axes = 0; axes < 3; axes ++){
          if(fieldtype == FIELD_D)      P[i].D[axes]       = ZA_D[axes][coord_q];
          if(fieldtype == FIELD_dDdy)   P[i].dDdy[axes]    = ZA_D[axes][coord_q];
          if(fieldtype == FIELD_ddDddy) P[i].D[axes]       = ZA_D[axes][coord_q]; // We store ddDddy in D to save memory
          if(fieldtype == FIELD_deltaD) P[i].dDdy[axes]    = ZA_D[axes][coord_q]; // We store deltaD in dDdy to save memory
        }
      } else if(LPTorder == 2){
        for(int axes = 0; axes < 3; axes ++){
          if(fieldtype == FIELD_D)      P[i].D2[axes]      = ZA_D[axes][coord_q];
          if(fieldtype == FIELD_dDdy)   P[i].dD2dy[axes]   = ZA_D[axes][coord_q];
          if(fieldtype == FIELD_ddDddy) P[i].D2[axes]      = ZA_D[axes][coord_q]; // We store ddD2ddy in D2 to save memory
          if(fieldtype == FIELD_deltaD) P[i].dD2dy[axes]   = ZA_D[axes][coord_q]; // We store deltaD2 in dD2dy to save memory
        }
      }

    } else {

      // Add to request buffer
      request_buffer[ cpu_id ][ n_ThisTask_needs_from_cpu[ cpu_id ] ].coord_q = coord_q;
      request_buffer[ cpu_id ][ n_ThisTask_needs_from_cpu[ cpu_id ] ].pid     = i;
      n_ThisTask_needs_from_cpu[ cpu_id ]++;
    }
  }

  // Show some info about how much data we need to communicate
  if(comm_verbose)
    for(int i = 0; i < NTask; i++)
      printf("Task [%3i] will send [%3i] disp_data to Task [%3i]\n", ThisTask, n_ThisTask_needs_from_cpu[i], i);

  // Allocate communication array
  // CPU #i will have to send n_to_send_all[i * NTask + j] pieces of particle data to CPU #j
  int *n_to_send_all = my_malloc(sizeof(int) * NTask * NTask);
  MPI_Allgather(&n_ThisTask_needs_from_cpu[0], NTask, MPI_INT, &n_to_send_all[0], NTask, MPI_INT, MPI_COMM_WORLD);

  // Print out after communication to check that all is fine 
  if(ThisTask == 0 && comm_verbose){
    for(int j = 0; j < NTask; j++){
      for(int i = 0; i < NTask; i++){
        printf("Task [%3i] will send [%3i] disp_data to Task [%3i]\n", j, n_to_send_all[j * NTask + i], i);
      }
    }
  }

  // Allocate send and recieve array
  struct disp_comm_data **request_data_to_send_to      = my_malloc(sizeof(struct disp_comm_data *) * NTask);
  struct disp_comm_data **request_data_to_recieve_from = my_malloc(sizeof(struct disp_comm_data *) * NTask);
  for(int i = 0; i < NTask; i++){
    request_data_to_send_to[i]      = my_malloc(sizeof(struct disp_comm_data) * n_to_send_all[ThisTask * NTask + i]);
    request_data_to_recieve_from[i] = my_malloc(sizeof(struct disp_comm_data) * n_to_send_all[ThisTask + NTask * i]);

    if(comm_verbose && (ThisTask == 0 && i != 0) ){
      printf("CPU [%3i] send to      [%3i] n = [%3i]\n", ThisTask, i, n_to_send_all[ThisTask * NTask + i]);
      printf("CPU [%3i] recieve from [%3i] n = [%3i]\n", ThisTask, i, n_to_send_all[ThisTask + NTask * i]);
    }

    // Reset recieve buffers
    for(int j = 0; j < n_to_send_all[ThisTask + NTask * i]; j++){
      request_data_to_recieve_from[i][j].pid     = 0;
      request_data_to_recieve_from[i][j].coord_q = 0;
    }

    // Assign the data we are about to send
    for(int j = 0; j < n_to_send_all[ThisTask * NTask + i]; j++){
      request_data_to_send_to[i][j].pid     = request_buffer[i][j].pid;
      request_data_to_send_to[i][j].coord_q = request_buffer[i][j].coord_q;
    }

    allocated_bytes += sizeof(struct disp_comm_data) * n_to_send_all[ThisTask * NTask + i];
    allocated_bytes += sizeof(struct disp_comm_data) * n_to_send_all[ThisTask + NTask * i];
  }

  // Free up some memory
  for(int i = 0; i < NTask; i++)
    my_free(request_buffer[i]);
  my_free(request_buffer);

  // Communicate
  for(int i = 1; i < NTask; i++){
    unsigned int send_request_to  = mymod(ThisTask + i, NTask);
    unsigned int get_request_from = mymod(ThisTask - i, NTask);
    unsigned int n_to_send    = n_to_send_all[ThisTask * NTask + send_request_to];
    unsigned int n_to_recieve = n_to_send_all[ThisTask + NTask * get_request_from];

    if(ThisTask == 0 && comm_verbose){
      printf("Send #[%i] to   #[%i]   Nsend    = [%i]\n", ThisTask, send_request_to,  n_to_send);
      printf("Get  #[%i] from #[%i]   Nrecieve = [%i]\n", ThisTask, get_request_from, n_to_recieve);
      printf("#[%i] Getting [%i] from #[%i]\n", ThisTask, n_to_recieve, get_request_from);
    }

    // Send to right, recieve from left
    MPI_Sendrecv(request_data_to_send_to[send_request_to], n_to_send * sizeof(struct disp_comm_data), MPI_BYTE, send_request_to, 0, 
        request_data_to_recieve_from[get_request_from], n_to_recieve * sizeof(struct disp_comm_data), MPI_BYTE, get_request_from, 0, MPI_COMM_WORLD, &status);
  }

  // Allocate return data memory
  struct disp_comm_data_return **return_data_to_send_to = my_malloc(sizeof(struct disp_comm_data_return *) * NTask);
  struct disp_comm_data_return **return_data_to_recieve_from = my_malloc(sizeof(struct disp_comm_data_return *) * NTask);
  for(int i = 0; i < NTask; i++){
    return_data_to_send_to[i]      = my_malloc(sizeof(struct disp_comm_data_return) * n_to_send_all[ThisTask + NTask * i]);
    return_data_to_recieve_from[i] = my_malloc(sizeof(struct disp_comm_data_return) * n_to_send_all[ThisTask * NTask + i]);

    allocated_bytes += sizeof(struct disp_comm_data_return) * n_to_send_all[ThisTask * NTask + i];
    allocated_bytes += sizeof(struct disp_comm_data_return) * n_to_send_all[ThisTask + NTask * i];
  }

  // Loop over all particles and assign the return_data
  for(int i = 0; i < NTask; i++){
    for(int j = 0; j < n_to_send_all[ThisTask + NTask * i]; j++){

      // Fetch request id's
      unsigned int coord_q = request_data_to_recieve_from[i][j].coord_q;
      unsigned int pid     = request_data_to_recieve_from[i][j].pid;

      // Compile up return data
      return_data_to_send_to[i][j].pid  = pid;
      for(int axes = 0; axes < 3; axes++){
        return_data_to_send_to[i][j].D[axes] = ZA_D[axes][coord_q];
      }
    }
  }

  // Communicate
  for(int i = 1; i < NTask; i++){
    unsigned int send_request_to  = mymod(ThisTask + i, NTask);
    unsigned int get_request_from = mymod(ThisTask - i, NTask);
    unsigned int n_to_send    = n_to_send_all[ThisTask + NTask * send_request_to];
    unsigned int n_to_recieve = n_to_send_all[ThisTask * NTask + get_request_from];

    if(ThisTask == 0 && comm_verbose){
      printf("Send #[%i] -> #[%i]   Nsend    = [%i]\n", ThisTask, send_request_to,  n_to_send);
      printf("Get  #[%i] <- #[%i]   Nrecieve = [%i]\n", ThisTask, get_request_from, n_to_recieve);
      printf("#[%i] Getting [%i] from #[%i]\n", ThisTask, n_to_recieve, get_request_from);
    }

    // Send to right, recieve from left
    MPI_Sendrecv(return_data_to_send_to[send_request_to], n_to_send * sizeof(struct disp_comm_data_return), MPI_BYTE, send_request_to, 0, 
        return_data_to_recieve_from[get_request_from], n_to_recieve * sizeof(struct disp_comm_data_return), MPI_BYTE, get_request_from, 0, MPI_COMM_WORLD, &status);
  }

  // Assign the displacement field to the remaining particles
  for(int i = 0; i < NTask; i++){
    for(int j = 0; j < n_to_send_all[ThisTask * NTask + i]; j++){
      unsigned int pid  = return_data_to_recieve_from[i][j].pid;
      float_kind *D = &(return_data_to_recieve_from[i][j].D[0]);

      // Check for out of bounds error just in case
      if(pid < 0 || pid > NumPart) printf("Error: %i %i\n", pid, NumPart);

      if(LPTorder == 1){
        for(int axes = 0; axes < 3; axes++){
          if(fieldtype == FIELD_D)      P[pid].D[axes]       = D[axes];
          if(fieldtype == FIELD_dDdy)   P[pid].dDdy[axes]    = D[axes];
          if(fieldtype == FIELD_ddDddy) P[pid].D[axes]       = D[axes]; // Save ddDddy in D to save memory
          if(fieldtype == FIELD_deltaD) P[pid].dDdy[axes]    = D[axes]; // Save deltaD in dDdy to save memory
        }
      } else if(LPTorder == 2){
        for(int axes = 0; axes < 3; axes++){
          if(fieldtype == FIELD_D)      P[pid].D2[axes]      = D[axes];
          if(fieldtype == FIELD_dDdy)   P[pid].dD2dy[axes]   = D[axes];
          if(fieldtype == FIELD_ddDddy) P[pid].D2[axes]      = D[axes]; // Save ddD2ddy in D2 to save memory
          if(fieldtype == FIELD_deltaD) P[pid].dD2dy[axes]   = D[axes]; // Save deltaD2 in dD2dy to save memory
        }
      }
    }
  }

  // Free up communication buffers
  for(int i = 0; i < NTask; i++){
    my_free(return_data_to_recieve_from[i] ); 
    my_free(return_data_to_send_to[i]      ); 
    my_free(request_data_to_recieve_from[i]);
    my_free(request_data_to_send_to[i]     );
  }
  my_free(return_data_to_recieve_from );
  my_free(return_data_to_send_to      );
  my_free(request_data_to_recieve_from);
  my_free(request_data_to_send_to     );
  my_free(n_ThisTask_needs_from_cpu   );
  my_free(n_to_send_all);

  // Free up displacment-fields
  for(int axes = 0; axes < 3; axes++) {
    my_free(ZA_D[axes]);
  }

  if(ThisTask == 0 && comm_verbose){
    printf("Memory used %f MB. For comparison density grid holds: %f MB\n", 
        allocated_bytes / 1024. / 1024., Total_size * sizeof(complex_kind) / 1024. / 1024.);
  }
}

#endif

