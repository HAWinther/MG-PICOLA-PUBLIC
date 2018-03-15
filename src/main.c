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

//=======================================================//
// This file contains the main driver routine for PICOLA //
//=======================================================//

#include "vars.h"
#include "proto.h"
#include "timer.h"
#include "msg.h"

int main(int argc, char **argv) {

  //======================================
  // Set up MPI
  //======================================
  ierr = MPI_Init(&argc, &argv); 
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  my_fftw_mpi_init();

  // Start timing
  msg_init();
  timer_set_category(_Init);

  //======================================
  // Read command line arguments
  //======================================
  if(argc < 2) {
    if(ThisTask == 0) {
      fprintf(stdout, "Input parameters not found\n");
      fprintf(stdout, "Call with <ParameterFile>\n");
    }
    ierr = MPI_Finalize();
    exit(1);
  }

  //======================================
  // Read the run parameters and setup code
  //======================================
  int i, stepDistr;
  if (ThisTask == 0) {
    printf("\n==============================================\n");
    printf("Reading Input Parameters and setting up PICOLA\n");
    printf("==============================================\n\n");
  }
  read_parameterfile(argv[1]);
  read_outputs();
  set_units();
#ifdef LIGHTCONE
  set_lightcone();
#endif

  if (UseCOLA){
    stepDistr   = 0;
    StdDA       = 0;
  } else{
    stepDistr   = 1;
    StdDA       = 2;
  }
  if (StdDA == 0){
    fullT = 1;
    nLPT  = -2.5;
  }

  //======================================
  // Show some info about the simulation
  //======================================
  if(ThisTask == 0) {
    printf("================\n");
    printf("Run Parameters\n");
    printf("================\n");
    printf("Cosmology:\n");
    printf("  Omega Matter(z=0) = %lf\n",Omega);
    printf("  Omega Baryon(z=0) = %lf\n",OmegaBaryon);
    printf("  Hubble Parameter(z=0) = %lf\n",HubbleParam);
    printf("  Sigma8(z=0) = %lf\n",Sigma8);
#ifndef GAUSSIAN
    printf("  F_nl = %lf\n",Fnl);
#endif
    printf("  Primordial Index = %lf\n",PrimordialIndex);
    printf("  Initial Redshift  = %lf\n",Init_Redshift);
#ifndef GAUSSIAN
    printf("  F_nl Redshift  = %lf\n",Fnl_Redshift);
#endif
    printf("\nSimulation:\n");
    printf("  Nmesh = %d\n", Nmesh);
    printf("  Nsample = %d\n", Nsample);
    printf("  Boxsize = %lf\n", Box);
    printf("  Buffer Size = %lf\n", Buffer);
    switch(WhichSpectrum) {
      case 0:
        switch (WhichTransfer) {
          case 1:
            printf("  Using Tabulated Transfer Function\n");
            break;
          default:
            printf("  Using Eisenstein & Hu Transfer Function\n");
            break;
        }
        break;
      case 1:
        printf("  Using Tabulated Power Spectrum\n");
        break;   
      default:
        printf("  Using Eisenstein & Hu Power Spectrum\n");
        break;
    }      
    if (UseCOLA) {
      printf("  Using COLA method\n");
    } else {
      printf("  Using Standard PM method\n");
    }
#ifdef LIGHTCONE
    printf("\nLightcone:\n");
    printf("  Maximum Comoving Radius = %lf\n", Light / Hubble * SphiStd(1.0 / (1.0+OutputList[0].Redshift), 1.0));
    printf("  Origin (x, y, z) = %lf, %lf, %lf\n", Origin_x, Origin_y, Origin_z);
    printf("  Nrep_min (x, y, z) = %d (%lf Mpc/h), %d (%lf Mpc/h), %d (%lf Mpc/h)\n", Nrep_neg_x, -Nrep_neg_x*Box-Origin_x, Nrep_neg_y, -Nrep_neg_y*Box-Origin_y, Nrep_neg_z, -Nrep_neg_z*Box-Origin_z);
    printf("  Nrep_max (x, y, z) = %d (%lf Mpc/h), %d (%lf Mpc/h), %d (%lf Mpc/h)\n", Nrep_pos_x, (Nrep_pos_x+1)*Box-Origin_x, Nrep_pos_y, (Nrep_pos_y+1)*Box-Origin_y, Nrep_pos_z, (Nrep_pos_z+1)*Box-Origin_z);
#endif
    printf("\nOutputs:\n");
    for (i = 0; i < Noutputs; i++) printf("  Redshift = %lf, Nsteps = %d\n", OutputList[i].Redshift, OutputList[i].Nsteps);
    fflush(stdout);
  }   

  if (ThisTask == 0) {
    printf("\n=============================================\n");
    printf(  "Initialising Transfer Function/Power Spectrum\n");
    printf(  "=============================================\n\n");
  }
  initialize_transferfunction();
  initialize_powerspectrum();
  initialize_ffts();
  initialize_parts();
  create_MPI_type_for_Particles(&PartDataMPIType);

  //=======================================================
  // Do the initialization of the modified gravity version
  // Compute growth-factors and make splines
  //=======================================================
  init_modified_version();

  timer_set_category(_GenerateIC);

  //===========================================================================================
  // Create the calculate the Zeldovich and 2LPT displacements and create the initial conditions
  //===========================================================================================
  int j, k, m;
  int NoutputStart = 0;
  int timeStep;
  unsigned int coord = 0;
  double da    = 0;
  double A     = 1.0/(1.0 + Init_Redshift);  // This is the scale factor which we'll be advancing below.
  double Di    = growth_D(A);                // initial growth factor
  double Di2   = growth_D2(A);               // initial 2nd order growth factor  
  double dDdy  = growth_dDdy(A);             // T[D_{za}]=dD_{za}/dy
  double dD2dy = growth_dD2dy(A);            // T[D_{2lpt}]=dD_{2lpt}/dy

  //===========================================================================================
  // A is the scale factor of the particle positions.
  // AI is the scale factor of the particle velocities.
  // AF is the scale factor to which we should kick the particle velocities.
  // AFF is the scale factor to which we should drift the particle positions.
  //===========================================================================================
  double AI = A, AF = A, AFF = A;  
  aexp_global = A;

  if(ReadParticlesFromFile){

    //===========================================================================================
    // Read particle from RAMSES, GADGET or ASCII file(s) and use it to create the displacement-field
    // All methods are in file [readICfromfile.c]
    // After density(k,zi) has been created we use LCDM growth-factor to bring it to z = 0 and
    // then call displacement_fields()
    //===========================================================================================

    P = my_malloc((size_t) (ceil(NumPart*Buffer)) * sizeof(struct part_data));

#ifdef READICFROMFILE
    ReadFilesMakeDisplacementField();
#endif

  } else {

    //===========================================================================================
    // Generate IC from scratch and Compute the displacementfield
    //===========================================================================================

    displacement_fields();
  }

  if(ThisTask == 0){
    printf("\n=================================\n");
    printf(  "Assign particle pos and vel...   \n");
    printf(  "=================================\n");
  }

  //===========================================================================================
  // Generate the initial particle positions and velocities
  // If UseCOLA = 0 (non-COLA), then velocity is ds/dy, which is simply the 2LPT IC.
  // Else set vel = 0 if we subtract LPT. This is the same as the action of the 
  // operator L_- from TZE, as initial velocities are in 2LPT.
  //===========================================================================================

  // Allocate memory for the particles
  if(!ReadParticlesFromFile){
    P = my_malloc((size_t) (ceil(NumPart*Buffer)) * sizeof(struct part_data));
  }

#ifdef SCALEDEPENDENT

  // Store the particle IDs 
  for(i = 0; i < Local_np; i++) {
    for (j = 0; j < Nsample; j++) {
      for (k = 0; k < Nsample; k++) {
        coord = (i * Nsample + j) * Nsample + k;
        P[coord].coord_q = coord;
        P[coord].init_cpu_id = ThisTask;
      }
    }
  }

  // Assign 1LPT and 2LPT displacementfields to the particles
  // Place [D] in Particles.D and [dDdy] in Particles.dDdy
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_D,      LPT_ORDER_ONE);
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_D,      LPT_ORDER_TWO);
  if(ThisTask == 0) printf("\n");

  assign_displacment_field_to_particles(A, AF, AFF, FIELD_dDdy,   LPT_ORDER_ONE);
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_dDdy,   LPT_ORDER_TWO);
  if(ThisTask == 0) printf("\n");

#endif

  // Assign positions, velocities and displacement
  for(i = 0; i < Local_np; i++) {
    for (j = 0; j < Nsample; j++) {
      for (k = 0; k < Nsample; k++) {
        coord = (i * Nsample + j) * Nsample + k;

#ifdef PARTICLE_ID          
        P[coord].ID = ((unsigned long long)((i + Local_p_start) * Nsample + j)) * (unsigned long long)Nsample + (unsigned long long)k;
#endif

        for (m = 0; m < 3; m++) {

#ifndef SCALEDEPENDENT
          //========================================
          // Assign displacementfields to particles
          //========================================
          P[coord].D[m]  = ZA[m][coord];
          P[coord].D2[m] = LPT[m][coord];
#endif

          if (UseCOLA == 0) {

            //==============================================================
            // Initial 2LPT velocity v = Psi^(1) * dD1/dt + Psi^(2) * dD2/dt
            //==============================================================  
#ifdef SCALEDEPENDENT
            P[coord].Vel[m] = P[coord].dDdy[m] + P[coord].dD2dy[m];
#else
            P[coord].Vel[m] = P[coord].D[m] * dDdy + P[coord].D2[m] * dD2dy;
#endif
          } else {

            P[coord].Vel[m] = 0.0;

          }
        }

        //======================================================================================================================
        // Initial 2LPT position: x = q + Psi^(1) D1 + Psi(2) D2 
        //======================================================================================================================

#ifdef SCALEDEPENDENT
        P[coord].Pos[0] = periodic_wrap((i + Local_p_start)*(Box / (double)Nsample) + P[coord].D[0] + P[coord].D2[0]);
        P[coord].Pos[1] = periodic_wrap( j                 *(Box / (double)Nsample) + P[coord].D[1] + P[coord].D2[1]);
        P[coord].Pos[2] = periodic_wrap( k                 *(Box / (double)Nsample) + P[coord].D[2] + P[coord].D2[2]);   
#else
        P[coord].Pos[0] = periodic_wrap((i + Local_p_start)*(Box / (double)Nsample) + P[coord].D[0] * Di + P[coord].D2[0] * Di2);
        P[coord].Pos[1] = periodic_wrap( j                 *(Box / (double)Nsample) + P[coord].D[1] * Di + P[coord].D2[1] * Di2);
        P[coord].Pos[2] = periodic_wrap( k                 *(Box / (double)Nsample) + P[coord].D[2] * Di + P[coord].D2[2] * Di2);   
#endif
        if(ThisTask == 0 && coord < 10) printf("Particle [%i] :   %8.3f   %8.3f   %8.3f\n", coord, P[coord].Pos[0], P[coord].Pos[1], P[coord].Pos[2]);
      }
    }
  }

  // Free up some memory
#ifndef SCALEDEPENDENT
  for (i = 0; i < 3; i++) {
    my_free(ZA[i]);
    my_free(LPT[i]);
  }
#endif

  // ===================================================================================================================
  // If we want to output or start the lightcone at the initial redshift this is where we do it (it is tricky to compare
  // floating point numbers due to rounding errors so instead we see whether they are close)
  // ===================================================================================================================
  if (((Init_Redshift-OutputList[0].Redshift)/Init_Redshift <= 1.0E-6) || (Init_Redshift <= 1.0e-6)) {

#ifndef LIGHTCONE

    if (ThisTask == 0) {
      printf("\n=============================\n");
      printf("Outputting Initial Conditions\n");
      printf("=============================\n\n");
    }

    for(int axes = 0; axes < 3; axes++)
      sumxyz[axes] = 0;

    // Output particles
    Output(A, AF, AFF, dDdy, dD2dy, OutputList[0].Redshift);

    // If this is the only output timestep then simply skip to the end
    if(Noutputs == 1) {
      goto finalize;
    }

#endif

    NoutputStart++;
  }

  //========================================================================================
  // Now, we get to the N-Body part where we evolve with time via the Kick-Drift-Kick Method
  //========================================================================================

  // The density grid and force grids and associated fftw plans
#ifndef MEMORY_MODE
  density = my_malloc(2 * Total_size * sizeof(float_kind));
  P3D     = (complex_kind *) density;

  N11     = my_malloc(2 * Total_size * sizeof(float_kind));
  N12     = my_malloc(2 * Total_size * sizeof(float_kind));
  N13     = my_malloc(2 * Total_size * sizeof(float_kind));
  FN11    = (complex_kind *) N11;
  FN12    = (complex_kind *) N12;
  FN13    = (complex_kind *) N13;

  plan    = my_fftw_mpi_plan_dft_r2c_3d(Nmesh, Nmesh, Nmesh, density, P3D, MPI_COMM_WORLD, FFTW_MEASURE); 
  p11     = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, FN11,    N11, MPI_COMM_WORLD, FFTW_MEASURE);
  p12     = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, FN12,    N12, MPI_COMM_WORLD, FFTW_MEASURE);
  p13     = my_fftw_mpi_plan_dft_c2r_3d(Nmesh, Nmesh, Nmesh, FN13,    N13, MPI_COMM_WORLD, FFTW_MEASURE);

  // Modified gravity allocation
  if(modified_gravity_active) AllocateMGArrays();
#endif

  if(ThisTask == 0) {
    printf("\n======================\n");
    printf("Beginning timestepping\n");
    printf("======================\n\n");
    fflush(stdout);
  }   
  timer_set_category(_TimeStepping);

  //=================================================
  // Loop over all the timesteps in the timestep list
  //=================================================
  timeSteptot = 0;

  // For light-cone sims we do not run the loop below for i = Noutputs
  int Nend = Noutputs;
#ifdef LIGHTCONE
  Nend = Noutputs-1;
#else
  Nend = Noutputs;
#endif
  for (i = NoutputStart; i <= Nend;i++){
    int nsteps = 0;
    double ao = 0;
    if (i == Noutputs) {
      nsteps = 1;
    } else {
      nsteps = OutputList[i].Nsteps;
      ao  = 1.0/(1.0+OutputList[i].Redshift);
      if (stepDistr == 0) da = (ao-A)/((double)nsteps);
      if (stepDistr == 1) da = (log(ao)-log(A))/((double)nsteps);
      if (stepDistr == 2) da = (CosmoTime(ao)-CosmoTime(A))/((double)nsteps);
    }

    //=========================================================
    // Perform the required number of timesteps between outputs
    //=========================================================
    for (timeStep = 0; timeStep < nsteps; timeStep++) {

      timeStep_global = timeStep;
      NoutputStart_global = NoutputStart;

#ifdef LIGHTCONE

      //=====================================================================================
      // For a lightcone simulation we always want the velocity set to mid-point of interval.
      //=====================================================================================

      if (stepDistr == 0) AF = A + da*0.5;
      if (stepDistr == 1) AF = A * exp(da*0.5);
      if (stepDistr == 2) AF = AofTime((CosmoTime(AFF)+CosmoTime(A))*0.5); 

#else
      //=====================================================================================
      // Calculate the time to update to
      // Half timestep for kick at output redshift to bring the velocities and positions to the same time
      //=====================================================================================

      if ((timeStep == 0) && (i != NoutputStart)) {
        AF = A; 
      } else { 

        //=====================================================================================
        // Set to mid-point of interval. In the infinitesimal timestep limit, these choices are identical. 
        // How one chooses the mid-point when not in that limit is really an extra degree of freedom in the code 
        // but Tassev et al. report negligible effects from the different choices below. 
        // Hence, this is not exported as an extra switch at this point.
        //=====================================================================================

        if (stepDistr == 0) AF = A + da*0.5;
        if (stepDistr == 1) AF = A * exp(da*0.5);
        if (stepDistr == 2) AF = AofTime((CosmoTime(AFF)+CosmoTime(A))*0.5); 
      }
#endif
      if (stepDistr == 0) AFF = A + da;
      if (stepDistr == 1) AFF = A * exp(da);
      if (stepDistr == 2) AFF = AofTime(CosmoTime(A)+da);

      // Compute the maximum particle number over all the CPUs
      int NumPart_max = NumPart;
      ierr = MPI_Allreduce(MPI_IN_PLACE, &NumPart_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      timeSteptot++;
      if (ThisTask == 0) {
        printf("Iteration = %d\n------------------\n",timeSteptot);
        printf("Maximum Particle Memory Used: %f %%\n", NumPart_max / (double)(Local_np*Nsample*Nsample*Buffer) * 100.0);
        if (i != Noutputs) {
          printf("a = %lf -> a = %lf\n", A, AFF);
          printf("z = %lf -> z = %lf\n", 1.0/A-1.0, 1.0/AFF-1.0);
          fflush(stdout);
        } else {
          printf("Final half timestep to update velocities...\n");
          fflush(stdout);
        }
      }

      print_memory_summary();

      //================================================================
      // Copy value of A to global variable. Needed by the mg.h routines
      //================================================================
      aexp_global = A;

      //=======================================================
      // Calculate the particle accelerations for this timestep
      //=======================================================
      GetDisplacements();

      /**********************************************************************************************/
      // If we wanted to interpolate the lightcone velocities we could put section currently in the // 
      // Drift subroutine here. This would allow us to have the velocities at AF and AFF which we   //
      // can use to invert the previous drift step and get the particle positions at AF and AFF     //
      /**********************************************************************************************/

#ifdef SCALEDEPENDENT
      // Compute the displacement fields for cases where we have scale depdendent growth
      if(ThisTask == 0){
        printf("\n=================================\n");
        printf("Compute displacment fields...   \n");
        printf("=================================\n");
      }

      // Place [ddDddy] in Particles.D and [deltaD] in Particles.dDdy
      assign_displacment_field_to_particles(A, AF, AFF, FIELD_deltaD, LPT_ORDER_ONE);
      assign_displacment_field_to_particles(A, AF, AFF, FIELD_deltaD, LPT_ORDER_TWO);
      if(ThisTask == 0) printf("\n");

      assign_displacment_field_to_particles(A, AF, AFF, FIELD_ddDddy, LPT_ORDER_ONE);
      assign_displacment_field_to_particles(A, AF, AFF, FIELD_ddDddy, LPT_ORDER_TWO);
      if(ThisTask == 0) printf("\n");

#endif
      //=============================
      // Kick the particle velocities
      //=============================
      if (ThisTask == 0) {
        printf("Kicking the particles...\n");
        fflush(stdout);
      }

      Kick(AI, AF, A, Di);

#ifndef LIGHTCONE

      //======================================================================================
      // If we are at an output timestep we modify the velocities and output the
      // particles. Then, if we are not yet at the end of the simulation, we update  
      // the velocity again up to the middle of the next timestep as per the usual KDK method.
      //======================================================================================

      if ((timeStep == 0) && (i != NoutputStart)) {

        if (ThisTask == 0) {
          printf("Outputting the particles...\n");
          fflush(stdout);
        }

        //======================================================================================
        // At the output timestep, add back LPT velocities if we had subtracted them. 
        // This corresponds to L_+ operator in TZE.
        //======================================================================================

        dDdy  = growth_dDdy(A);    // dD_{za}/dy
        dD2dy = growth_dD2dy(A);   // dD_{2lpt}/dy

        Output(A, AF, AFF, dDdy, dD2dy, OutputList[i-1].Redshift);

        //======================================================================================
        // If we have reached the last output timestep we skip to the end
        //======================================================================================

        if(i == Noutputs) {
          if (ThisTask == 0) {
            printf("Iteration %d finished\n------------------\n\n", timeSteptot);
            fflush(stdout);
          }

          // Remember to free up displacement-field before we exit
          for (j = 0; j < 3; j++) my_free(Disp[j]);

          goto finalize;
        }

        // =====================================================================================
        // Otherwise we simply update the velocity to the middle of the timestep, where it
        // would have been if we weren't outputting. This involves only recalculating and applying `dda'
        // as the acceleration remains the same as calculated earlier
        // =====================================================================================

        AI = A;
        if (stepDistr == 0) AF = A + da * 0.5;
        if (stepDistr == 1) AF = A * exp(da * 0.5);
        if (stepDistr == 2) AF = AofTime((CosmoTime(AFF)+CosmoTime(A)) * 0.5); 
        for(int axes = 0; axes < 3; axes++)
          sumDxyz[axes] = 0;

        Kick(AI, AF, A, Di);
      }
#endif

      // Free up displacement-field
      for (j = 0; j < 3; j++) my_free(Disp[j]);

      // ============================
      // Drift the particle positions
      // ============================

      if (ThisTask == 0) {
        printf("Drifting the particles...\n");
        fflush(stdout);
      }

#ifdef LIGHTCONE
      if (i > 0) {
        Drift_Lightcone(A, AFF, AF, Di, Di2);
      } else {
        Drift(A, AFF, AF, Di, Di2);
      }
#else
      Drift(A, AFF, AF, Di, Di2);
#endif
      //=================
      // Step in time
      //=================

      A  = AFF;
      AI = AF;

      Di  = growth_D(A);
      Di2 = growth_D2(A);

      if (ThisTask == 0) {
        printf("Iteration %d finished\n------------------\n\n", timeSteptot);
        fflush(stdout);
      }

      ierr = MPI_Barrier(MPI_COMM_WORLD); 
    }
  }

  //============================
  // Here is the last little bit
  //============================
#ifndef LIGHTCONE
finalize:
#endif

  if (ThisTask == 0) {
    printf("===============================\n");
    printf("Finishing up! Timing:\n");
    printf("===============================\n");
    fflush(stdout);
  }

#ifdef LIGHTCONE
  Output_Info_Lightcone();
#endif

  // Free up memory
  free_powertable();
  free_transfertable();

  my_free(P);
  my_free(OutputList);
  my_free(Slab_to_task);
  my_free(Part_to_task);
  my_free(Local_nx_table);
  my_free(Local_np_table);

#ifdef GENERIC_FNL
  my_free(KernelTable);
#endif

#ifdef LIGHTCONE
  my_free(Noutput);
  my_free(repflag);
#endif

#ifndef MEMORY_MODE
  my_free(density);
  my_free(N11);
  my_free(N12);
  my_free(N13);  
  my_fftw_destroy_plan(plan);
  my_fftw_destroy_plan(p11);
  my_fftw_destroy_plan(p12);
  my_fftw_destroy_plan(p13);
  if(modified_gravity_active) FreeMGArrays();
#endif

  free_up_splines();

#ifdef MASSIVE_NEUTRINOS
  free_CAMB_splines();
#endif

#ifdef SCALEDEPENDENT
  my_free(cdelta_cdm);
  my_free(cdelta_cdm2);
#endif

  my_fftw_mpi_cleanup();

  timer_print();

  print_memory_summary();

  MPI_Finalize();    

  return 0;
}

//================================
// Kicking the particle velocities
//================================
void Kick(double AI, double AF, double A, double Di) {
  timer_start(_Kick);
  double dda;
  double force[3];

  for(int axes = 0; axes < 3; axes++)
    sumxyz[axes] = 0;

  if (StdDA == 0) {
    dda = Sphi(AI, AF, A);
  } else if (StdDA == 1) {
    dda = (AF - AI) * A / Qfactor(A);
  } else {
    dda = SphiStd(AI, AF);
  }  

#ifdef SCALEDEPENDENT

  // Update velocity
  for(unsigned int n = 0; n < NumPart; n++) {
    for(int axes = 0; axes < 3; axes ++) {
      Disp[axes][n]  -= sumDxyz[axes];

      // NB: D/D2 here contains ddDddy/ddD2ddy. This is done to save memory
      force[axes]     = -1.5 * Omega * Disp[axes][n] - UseCOLA * ( P[n].D[axes] + P[n].D2[axes] ) / A;
      P[n].Vel[axes] += force[axes] * dda;
      sumxyz[axes]   += P[n].Vel[axes];
    }
  }

#else

  // Second derivate of the growth factors
  double ddDddy  = growth_ddDddy(A);   // T^2[D_{ZA}]=d^2 D_{ZA}/dy^2
  double ddD2ddy = growth_ddD2ddy(A);  // T^2[D_{2lpt}]=d^2 D_{2lpt}/dy^2

  // Update velocity
  for(unsigned int n = 0; n < NumPart; n++) {
    for(int axes = 0; axes < 3; axes ++) {
      Disp[axes][n]  -= sumDxyz[axes];
      force[axes]     = -1.5 * Omega * Disp[axes][n] - UseCOLA * ( P[n].D[axes] * ddDddy + P[n].D2[axes] * ddD2ddy ) / A;
      P[n].Vel[axes] += force[axes] * dda;
      sumxyz[axes]   += P[n].Vel[axes];
    }
  }
#endif

  // Make sumx, sumy and sumz global averages
  for(int axes = 0; axes < 3; axes ++){
    ierr = MPI_Allreduce(MPI_IN_PLACE, &(sumxyz[axes]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    sumxyz[axes] /= (double) TotNumPart;
  }

  timer_stop(_Kick);
}

//================================
// Drifting the particle positions
//================================
void Drift(double A, double AFF, double AF, double Di, double Di2) {
  timer_start(_Drift);
  double dyyy;

  if (StdDA == 0) {
    dyyy = Sq(A, AFF, AF);
  } else if (StdDA == 1) {
    dyyy = (AFF - A) / Qfactor(AF);
  } else {
    dyyy = SqStd(A, AFF);
  }

#ifdef SCALEDEPENDENT

  // Update positions
  for(unsigned int n = 0; n < NumPart; n++) {
    for(int axes = 0; axes < 3; axes++){
      P[n].Pos[axes] += (P[n].Vel[axes] - sumxyz[axes]) * dyyy;

      // NB: dDdy here contain (deltaD = D(AFF) - D(A)). This is done to save memory
      P[n].Pos[axes]  = periodic_wrap(P[n].Pos[axes] + UseCOLA*( P[n].dDdy[axes] + P[n].dD2dy[axes] ));
    }
  }

#else

  // Change in growth-factors
  double deltaD  = (growth_D(AFF) - Di);    // change in D
  double deltaD2 = (growth_D2(AFF) - Di2);  // change in D_{2lpt}

  // Update positions
  for(unsigned int n = 0; n < NumPart; n++) {
    for(int axes = 0; axes < 3; axes++){
      P[n].Pos[axes] += (P[n].Vel[axes] - sumxyz[axes]) * dyyy;
      P[n].Pos[axes]  = periodic_wrap(P[n].Pos[axes] + UseCOLA*(P[n].D[axes] * deltaD + P[n].D2[axes] * deltaD2 ));
    }
  }
#endif

  timer_stop(_Drift);
}

//=================
// Output the data
//=================
void Output(double A, double AF, double AFF, double dDdy, double dD2dy, double z_from_outputlist) {
  timer_set_category(_Output);
  FILE * fp; 
  char buf[300];
  int nprocgroup, groupTask, masterTask;
  unsigned int n;
  double Z         = z_from_outputlist;
  double fac       = Hubble / pow(A,1.5);
  double lengthfac = UnitLength_in_cm / 3.085678e24;     // Convert positions to Mpc/h
  double velfac    = UnitVelocity_in_cm_per_s / 1.0e5;   // Convert velocities to km/s

  // Are we outputting at the first step?
  int first_step = ( fabs(Init_Redshift - Z)/Init_Redshift <= 1.0E-6 ); 
  (void) first_step;

#ifdef SCALEDEPENDENT

  /* We can avoid some computation below at the expense of copying
   * but this requires memory
   * We make a copy of the displacment-field [dDdy] which contains [deltaD]
   * and add in dDdy in this field needed for velocities below
   float *tmp_buffer_deltaD  = my_malloc(sizeof(float) * NumPart * 3);
   float *tmp_buffer_deltaD2 = my_malloc(sizeof(float) * NumPart * 3);
   for(unsigned int i = 0; i < NumPart; i++){
   for(int axes = 0; axes < 3; axes++){
   tmp_buffer_deltaD[3*i + axes]  = P[i].dDdy[axes];
   tmp_buffer_deltaD2[3*i + axes] = P[i].dD2dy[axes];
   }
   }
   */

  // Compute dDdy
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_dDdy, LPT_ORDER_ONE);
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_dDdy, LPT_ORDER_TWO);
  if(ThisTask == 0) printf("\n");

#endif

#ifdef MATCHMAKER_HALOFINDER
  //==================================================================================================================
  // Run FoF halofinder. NB: For scaledependent case we must have dD/dy and dD2/dy in the P->dDdy and P->dD2dy vectors
  // This is the case here. As currently written we duplicate particles in this routine which requires extra memory.
  //==================================================================================================================
  if(mm_run_matchmaker && !first_step) {
    struct PicolaToMatchMakerData data;
    data.output_format  = mm_output_format;
    data.output_pernode = mm_output_pernode;
    data.np_min         = mm_min_npart_halo;
    data.b_fof          = mm_linking_length;
    data.NumPart        = NumPart;
    data.n_part_1d      = Nsample;
    data.omega_m        = Omega;
    data.omega_l        = 1.0 - Omega;
    data.HubbleParam    = HubbleParam;
    data.redshift       = 1.0/A - 1.0;
    data.Local_p_start  = (int)Local_p_start;
    data.P              = P;
    data.growth_dDdy    = &growth_dDdy;
    data.growth_dD2dy   = &growth_dD2dy;

    // Particle mass in 10^{10} Msun/h (10^{10} = MATCHMAKER_MASS_FACTOR)
    data.mass_part      = (3.0 * Omega * Hubble * Hubble * Box * Box * Box) / (8.0 * PI * G * TotNumPart);

    // We use user-units in MatchMaker
    data.norm_pos       = lengthfac;
    data.boxsize        = Box * lengthfac;
    data.dx_extra       = mm_dx_extra_mpc * lengthfac;

    // We work with comoving velocity [ dx/dt ] in km/s in MatchMaker
    data.norm_vel       = (Hubble / A / A);

    sprintf(data.OutputDir, OutputDir);
    sprintf(data.FileBase,  FileBase);

    // Could perhaps try to free some memory here just to make sure we don't run out? 
    // If not in memory mode then we can do this. Consider adding this here

    timer_start(_HaloFinding);
    MatchMaker(data);
    timer_stop(_HaloFinding);
  }
#endif

#ifdef COMPUTE_POFK
  // Compute the RSD power-spectrum
  if(pofk_compute_rsd_pofk >= 1 && !first_step){
    compute_RSD_powerspectrum(A, 1);
  }
#endif

  timer_start(_WriteOutput);

#ifdef GADGET_STYLE
  size_t bytes;
  int k, pc, dummy, blockmaxlen;
  float * block;
#ifdef PARTICLE_ID
  unsigned long long * blockid;
#endif
#endif

  nprocgroup = NTask / NumFilesWrittenInParallel;
  if (NTask % NumFilesWrittenInParallel) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if(NumPart > 0) {
        int Zint, Zfrac;
        if(Z < 0.0){
          Zint = Zfrac = 0;
        } else {
          Zint  = (int)floor(Z);
          Zfrac = (int)((Z - Zint)*1000);
        }

        sprintf(buf, "%s/%s_z%dp%03d.%d", OutputDir, FileBase, Zint, Zfrac, ThisTask);

        if(!(fp = fopen(buf, "w"))) {
          printf("\nERROR: Can't write in file '%s'.\n\n", buf);
          FatalError((char *)"main.c Output, can't write to file!");
        }
        fflush(stdout);

#ifdef GADGET_STYLE

        //===============================================================
        // Write a GADGET file
        //===============================================================

        //===================
        // Write header
        //===================
        write_gadget_header(fp, A);

#ifdef MEMORY_MODE
        // Allocate buffer for output. We make it so that we have to gather and write 6 times
        bytes = sizeof(float) * NumPart;
        block = my_malloc(bytes);
#else
        // When not in memory mode we are a bit more conservative here
        bytes = 10*1024*1024;
        block = my_malloc(bytes);
#endif
        blockmaxlen = bytes / (3 * sizeof(float));

        //===================
        // Write coordinates
        //===================
        dummy = sizeof(float) * 3 * NumPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          for(k = 0; k < 3; k++) block[3 * pc + k] = (float)(lengthfac*P[n].Pos[k]);
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3 * pc, fp);
            pc = 0;
          }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3 * pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        //===================
        // Write velocities
        //===================
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          // Remember to add the ZA and 2LPT velocities back on and convert to PTHalos velocity units
#ifdef SCALEDEPENDENT
          for(k = 0; k < 3; k++) 
            block[3 * pc + k] = (float)(velfac*fac*(P[n].Vel[k] - sumxyz[k] 
                  + (P[n].dDdy[k] + P[n].dD2dy[k] ) * UseCOLA));
#else
          for(k = 0; k < 3; k++) 
            block[3 * pc + k] = (float)(velfac*fac*(P[n].Vel[k] - sumxyz[k] 
                  + (P[n].D[k] * dDdy + P[n].D2[k] * dD2dy ) * UseCOLA));
#endif
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3 * pc, fp);
            pc = 0;
          }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3 * pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

#ifdef PARTICLE_ID
        blockid = (unsigned long long *)block;
        blockmaxlen = bytes / sizeof(unsigned long long);

        //===================
        // Write particle ID
        //===================
        dummy = sizeof(unsigned long long) * NumPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(n = 0, pc = 0; n < NumPart; n++) {
          blockid[pc] = P[n].ID;
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(blockid, sizeof(unsigned long long), pc, fp);
            pc = 0;
          }
        }
        if(pc > 0) my_fwrite(blockid, sizeof(unsigned long long), pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
#endif
        my_free(block);   
#else
        //==========================================================
        // Output as ASCII
        //==========================================================

        // Total number of particle in each file on the first line
        fprintf(fp,"%u\n", NumPart);

        for(n = 0; n < NumPart; n++){
          double P_Vel[3];
          for(int axes = 0; axes < 3; axes++) {
#ifdef SCALEDEPENDENT
            P_Vel[axes] = fac*(P[n].Vel[axes] - sumxyz[axes] + (P[n].dDdy[axes] + P[n].dD2dy[axes] ) * UseCOLA);
#else
            P_Vel[axes] = fac*(P[n].Vel[axes] - sumxyz[axes] + (P[n].D[axes] * dDdy + P[n].D2[axes] * dD2dy) * UseCOLA);
#endif
          }

          // Output positions in Mpc/h and velocities in km/s
#ifdef PARTICLE_ID
          fprintf(fp,"%12llu %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
              P[n].ID, (float)(lengthfac*P[n].Pos[0]), (float)(lengthfac*P[n].Pos[1]), (float)(lengthfac*P[n].Pos[2]),
              (float)(velfac*P_Vel[0]),      (float)(velfac*P_Vel[1]),       (float)(velfac*P_Vel[2]));
#else
          fprintf(fp,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
              (float)(lengthfac*P[n].Pos[0]), (float)(lengthfac*P[n].Pos[1]), (float)(lengthfac*P[n].Pos[2]),
              (float)(velfac*P_Vel[0]),       (float)(velfac*P_Vel[1]),       (float)(velfac*P_Vel[2]));
#endif
        }
#endif
        fclose(fp);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD); 
  }
  Output_Info(A);

#ifdef SCALEDEPENDENT

  /* We can avoid some computation below at the expense of copying
   * but this requires memory
   * We make a copy of the displacment-field [dDdy] which contains [deltaD]
   * and add in dDdy in this field needed for velocities below
   for(unsigned int i = 0; i < NumPart; i++){
   for(int axes = 0; axes < 3; axes++){
   P[i].dDdy[axes]  = tmp_buffer_deltaD[3*i + axes];
   P[i].dD2dy[axes] = tmp_buffer_deltaD2[3*i + axes];
   }
   }
   my_free(tmp_buffer_deltaD);
   my_free(tmp_buffer_deltaD2);
   */

  // Compute deltaD
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_deltaD, LPT_ORDER_ONE);
  assign_displacment_field_to_particles(A, AF, AFF, FIELD_deltaD, LPT_ORDER_TWO);
  if(ThisTask == 0) printf("\n");

#endif

  timer_stop(_WriteOutput);
  timer_set_category(_TimeStepping);
  return;
}

//=====================================================================================
// Generate the info file which contains a list of all the output files, 
// the 8 corners of the slices on those files and the number of particles in the slice
//=====================================================================================
void Output_Info(double A) {
  FILE * fp; 
  char buf[300];
  double Z = (1.0/A) - 1.0;

  int *Local_p_start_table    = my_malloc(sizeof(int) * NTask);
  unsigned int *Noutput_table = my_malloc(sizeof(unsigned int) * NTask);
  int local_p_start_int = (int)Local_p_start;
  MPI_Allgather(&local_p_start_int, 1, MPI_INT, Local_p_start_table, 1, MPI_INT, MPI_COMM_WORLD); 
  MPI_Allgather(&NumPart, 1, MPI_UNSIGNED, Noutput_table, 1, MPI_UNSIGNED, MPI_COMM_WORLD); 

  if (ThisTask == 0) { // Added info_ to filename
    sprintf(buf, "%s/info_%s_z%dp%03d.txt", OutputDir, FileBase, (int)Z, (int)rint((Z-(int)Z)*1000));
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"main.c Output_Info, can't write to file");
    }
    fflush(stdout);
    fprintf(fp, "#    FILENUM      XMIN         YMIN        ZMIN         XMAX         YMAX         ZMAX         NPART    \n");
    double y0 = 0.0, z0 = 0.0, y1 = Box, z1 = Box;
    for (int i = 0; i < NTask; i++) {
      double x0 = Local_p_start_table[i] * (Box / (double)Nsample); 
      double x1 = (Local_p_start_table[i] + Local_np_table[i]) * (Box / (double)Nsample);
      fprintf(fp, "%12d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12u\n", i, x0, y0, z0, x1, y1, z1, Noutput_table[i]);
    }
    fclose(fp);
  }

  my_free(Noutput_table);
  my_free(Local_p_start_table);

  return;
}

//====================================================
// Write a GADGET header
//====================================================

#ifdef GADGET_STYLE
void write_gadget_header(FILE *fp, double A){
  double Z = 1.0/A - 1.0;
  int dummy;      

  // Gadget header stuff
  for(int k = 0; k < 6; k++) {
    header.npart[k]      = 0;
    header.npartTotal[k] = 0;
    header.mass[k]       = 0;
  }
  
  header.npart[1]      = NumPart;
  header.npartTotal[1] = (unsigned int) TotNumPart;
  header.npartTotalHighWord[1] = (unsigned int) (TotNumPart >> 32);
  header.mass[1]       = (3.0 * Omega * Hubble * Hubble * Box * Box * Box) / (8.0 * PI * G * TotNumPart);

  header.time          = A;
  header.redshift      = Z;

  header.flag_sfr        = 0;
  header.flag_feedback   = 0;
  header.flag_cooling    = 0;
  header.flag_stellarage = 0;
  header.flag_metals     = 0;
  header.flag_stellarage = 0;
  header.flag_metals     = 0;
  header.flag_entropy_instead_u = 0;

  header.num_files = NTaskWithN;

  // We might rescale the positions when outputting
  // so do the same to the Boxsize
  double lengthfac    = UnitLength_in_cm / 3.085678e24;
  header.BoxSize      = Box * lengthfac;
  
  header.Omega0       = Omega;
  header.OmegaLambda  = 1.0 - Omega;
  header.HubbleParam  = HubbleParam;

  dummy = sizeof(header);
  my_fwrite(&dummy,  sizeof(dummy),  1, fp);
  my_fwrite(&header, sizeof(header), 1, fp);
  my_fwrite(&dummy,  sizeof(dummy),  1, fp);
}
#endif
