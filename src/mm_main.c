///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of MatchMaker.                                  //
// This is a hacked version by Hans Winther                          //
//                                                                   //
// MatchMaker is free software: you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License as    //
// published by the Free Software Foundation, either version 3 of    //
// the License, or (at your option) any later version.               //
//                                                                   //
// MatchMaker is distributed in the hope that it will be useful, but //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with MatchMaker.  If not, see                               //
//        <http://www.gnu.org/licenses/>.                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include "mm_common.h"

//===================================================
// PICOLA header. We need particle struct definition
//===================================================
#include "vars.h"      

// We only need to MPI initialize MatchMaker once
int mm_mpi_is_initialized = 0;

static int mm_mpi_init(){
  //Define MPI Particle and Halo structures
  MPI_Aint p_off[5];
  MPI_Datatype p_ot[5];
  int p_bc[5];
  
  p_bc[0] = 3; //3*x
  p_bc[1] = 3; //3*v
  p_bc[2] = 1; //id
  p_bc[3] = 1; //fof_id
  p_bc[4] = 1; //cll_id
  
  p_off[0] = offsetof(Particle,x);
  p_off[1] = offsetof(Particle,v);
  p_off[2] = offsetof(Particle,id);
  p_off[3] = offsetof(Particle,fof_id);
  p_off[4] = offsetof(Particle,cll_id);
  
  p_ot[0] = MPI_FLOAT;
  p_ot[1] = MPI_FLOAT;
#ifdef MATCHMAKER_LONGIDS
  p_ot[2] = MPI_UNSIGNED_LONG_LONG;
#else //MATCHMAKER_LONGIDS
  p_ot[2] = MPI_UNSIGNED;
#endif //MATCHMAKER_LONGIDS
#ifdef MATCHMAKER_LONGINT
  p_ot[3] = MPI_LONG;
  p_ot[4] = MPI_LONG;
#else //MATCHMAKER_LONGINT
  p_ot[3] = MPI_INT;
  p_ot[4] = MPI_INT;
#endif //MATCHMAKER_LONGINT
  
  MPI_Type_struct(5,p_bc,p_off,p_ot,&ParticleMPI);
  
  MPI_Type_commit(&ParticleMPI);

  MPI_Aint h_off[12];
  MPI_Datatype h_ot[12];
  int h_bc[12];
  
  h_bc[0]  = 1; // np
  h_bc[1]  = 1; // m
  h_bc[2]  = 3; // 3*x_avg
  h_bc[3]  = 3; // 3*x_rms
  h_bc[4]  = 3; // 3*v_avg
  h_bc[5]  = 3; // 3*v_rms
  h_bc[6]  = 3; // 3*lam
  h_bc[7]  = 1; // b
  h_bc[8]  = 1; // c
  h_bc[9]  = 3; // 3*ea
  h_bc[10] = 3; // 3*ea
  h_bc[11] = 3; // 3*ea
  
  h_off[0]  = offsetof(FoFHalo,np);
  h_off[1]  = offsetof(FoFHalo,m_halo);
  h_off[2]  = offsetof(FoFHalo,x_avg);
  h_off[3]  = offsetof(FoFHalo,x_rms);
  h_off[4]  = offsetof(FoFHalo,v_avg);
  h_off[5]  = offsetof(FoFHalo,v_rms);
  h_off[6]  = offsetof(FoFHalo,lam);
  h_off[7]  = offsetof(FoFHalo,b);
  h_off[8]  = offsetof(FoFHalo,c);
  h_off[9]  = offsetof(FoFHalo,ea);
  h_off[10] = offsetof(FoFHalo,eb);
  h_off[11] = offsetof(FoFHalo,ec);
  
  h_ot[0]  = MPI_INT;
  h_ot[1]  = MPI_FLOAT;
  h_ot[2]  = MPI_FLOAT;
  h_ot[3]  = MPI_FLOAT;
  h_ot[4]  = MPI_FLOAT;
  h_ot[5]  = MPI_FLOAT;
  h_ot[6]  = MPI_FLOAT;
  h_ot[7]  = MPI_FLOAT;
  h_ot[8]  = MPI_FLOAT;
  h_ot[9]  = MPI_FLOAT;
  h_ot[10] = MPI_FLOAT;
  h_ot[11] = MPI_FLOAT;
  
  MPI_Type_struct(12,h_bc,h_off,h_ot,&HaloMPI);
  MPI_Type_commit(&HaloMPI);

  return 0;
}

void MatchMaker(struct PicolaToMatchMakerData data){
  // Initialize the MPI part of MatchMaker
  // Only needed the first time we call this routine
  if(mm_mpi_is_initialized == 0){
    mm_mpi_init();
    mm_msg_init();
    mm_mpi_is_initialized = 1;
  }

  mm_msg_printf("\n=========================\n");
  mm_msg_printf("Powering up MatchMaker   \n");
  mm_msg_printf("=========================\n\n");

  //============================================
  // Get parameters from PICOLA 
  //============================================
  Param.output_format   = data.output_format;
  Param.output_pernode  = data.output_pernode;
  Param.dx_extra        = data.dx_extra;
  Param.np_min          = data.np_min;
  Param.b_fof           = data.b_fof;
  Param.n_part_1d       = data.n_part_1d;
  Param.n_part          = (lint)data.n_part_1d * (lint)data.n_part_1d * (lint)data.n_part_1d;
  Param.boxsize         = data.boxsize;
  Param.omega_m         = data.omega_m;
  Param.omega_l         = data.omega_l;
  Param.redshift        = data.redshift;
  Param.norm_vel        = data.norm_vel;
  Param.norm_pos        = data.norm_pos;
  Param.h               = data.HubbleParam;
  Param.NumPart         = NumPart;
  Param.Local_p_start   = data.Local_p_start;
  Param.P               = data.P;
  Param.growth_dDdy     = data.growth_dDdy;
  Param.growth_dD2dy    = data.growth_dD2dy;
  sprintf(Param.OutputDir,"%s", data.OutputDir);
  sprintf(Param.FileBase, "%s", data.FileBase);
  Param.mp              = data.mass_part;
 
  // Check that the particle mass is OK
  double mass_expected = 27.7549954*Param.omega_m*pow(Param.boxsize,3)/Param.n_part;
  if( fabs(mass_expected / Param.mp - 1.0) > 0.01){
    mm_msg_printf("Warning: particle mass differs from expected %f != %f\n", mass_expected, Param.mp);
  }

  // Make output name
  double z  = Param.redshift;
  int zint  = (int)z;
  int zfrac = (int)((z - zint)*1000);
  sprintf(Param.output_prefix, "%s/matchmaker_%s_z%d.%03d", Param.OutputDir, Param.FileBase, zint, zfrac);

  // Get MPI parameters
  MPI_Comm_rank(MPI_COMM_WORLD,&(Param.i_node));
  MPI_Comm_size(MPI_COMM_WORLD,&(Param.n_nodes));
  Param.i_node_left = (Param.i_node-1+Param.n_nodes)%Param.n_nodes;
  Param.i_node_right = (Param.i_node+1)%Param.n_nodes;
  
  // Check buffer size
  if(Param.dx_extra >= 0.5 * Param.boxsize/(double) Param.n_nodes){
    mm_msg_printf("Buffer size might be too big! MatchMaker will not be run!\nReduce it or the number of cores!\n");
    mm_msg_printf("The maximum number of cores for this buffersize is: %i\n",  (int)(0.5 * Param.boxsize/Param.dx_extra));
    return;
  }

  //===============================================================
  // This is critical. If the number of tasks does not divide Nmesh
  // and/or Nsample then MatchMaker uses a different split of the
  // domain as PICOLA. We solve this by communicating from PICOLA
  // where the domain starts/ends for each task. Apart from 
  // x_offset we also need to recompute dx_domain when getting
  // particles from the right
  //===============================================================

  // Compute offset (position of left boundary) from PICOLA p-table
  Param.x_offset = Local_p_start * (Param.boxsize / (double)Param.n_part_1d); 

  // Get offset from the right and send to the left
  int tag = 100;
  MPI_Status stat;
  MPI_Sendrecv(&Param.x_offset,1,MPI_FLOAT,Param.i_node_left ,tag,
      &Param.x_offset_right,1,MPI_FLOAT,Param.i_node_right,tag,
      MPI_COMM_WORLD,&stat);

  // Compute size of current domain dx_domain 
  Param.dx_domain = Param.x_offset_right - Param.x_offset;
  if(Param.dx_domain < 0.0) Param.dx_domain += Param.boxsize;

#ifdef MATCHMAKER_DEBUG
  printf("Task [%i]  Holds positions %f < x < %f  dx_domain = %f\n", Param.i_node, Param.x_offset, Param.x_offset_right, Param.dx_domain);
#endif

  mm_msg_printf("Allocating particles and translating them to MatchMaker\n");
  Particles *particles = my_malloc(sizeof(Particles));
  picola_to_matchmaker_particles(particles);

  mm_msg_printf("Getting halos\n");
  lint n_halos;
  FoFHalo *fh = fof_get_halos(&n_halos,particles);

  mm_msg_printf("Writing output\n");
  write_halos(n_halos,fh);
  mm_msg_printf("=========================\n\n");

  // Clean up all memory
  my_free(particles->p);
  my_free(particles->p_back);
  my_free(particles);
}

//============================================================
// Compare particles according to their x coordinate
//============================================================
int compare_parts(const void *e1,const void *e2){
  const Particle *p1=e1;
  const Particle *p2=e2;

  if(p1->x[0]<p2->x[0])
    return -1;
  else if(p1->x[0]>p2->x[0])
    return 1;
  else
    return 0;
}

//============================================================
// Translating PICOLA particles to MatchMaker particles
// We count how many particles we will send/recieve and
// only allocate as much memory as we really need
//============================================================
void picola_to_matchmaker_particles(Particles *particles){
  lint i, j;
  int tag = 100;
  MPI_Status stat;
 
  // Count how many particles to send to the left
  lint np_toleft = 0;
  float edge_total_left = Param.x_offset;
  for(j=0;j<Param.NumPart;j++) {
    float x = Param.norm_pos * Param.P[j].Pos[0] - edge_total_left;
    if(x <= Param.dx_extra){
      np_toleft++;
    }
  }
  
  // Get how many to recieve from the right
  lint np_fromright;
#ifdef MATCHMAKER_LONGINT
  MPI_Sendrecv(&np_toleft,   1,MPI_LONG,Param.i_node_left ,tag,
      &np_fromright,1,MPI_LONG,Param.i_node_right,tag,
      MPI_COMM_WORLD,&stat);
#else //MATCHMAKER_LONGINT
  MPI_Sendrecv(&np_toleft,   1,MPI_INT,Param.i_node_left ,tag,
      &np_fromright,1,MPI_INT,Param.i_node_right,tag,
      MPI_COMM_WORLD,&stat);
#endif //MATCHMAKER_LONGINT
  
  // We now know how many particles we will have in the current node
  particles->np_fromright = np_fromright;
  particles->np_toleft    = np_toleft;
  particles->np_indomain  = Param.NumPart;
  particles->np_local     = particles->np_indomain+np_fromright;
  particles->np_allocated = particles->np_local;
  particles->np_centre    = particles->np_indomain-np_toleft;
  particles->np_total     = Param.n_part;
  particles->np_average   = (float) Param.n_part / (float) Param.n_nodes;
  particles->p_back       = my_malloc(particles->np_toleft*sizeof(Particle));
  particles->p            = my_malloc(sizeof(Particle)*particles->np_local);

  //==============================
  // Allocate partices
  //==============================
  unsigned long long np_tot = Param.n_part;
  float dx_domain           = Param.dx_domain;
  int this_node             = Param.i_node;

  // Translate from code vel to true velocity
  float norm_vel = Param.norm_vel;

  // Translate from code pos to pos in Mpc/h
  float norm_pos = Param.norm_pos;

#ifndef SCALEDEPENDENT
  double A     = 1.0 / (1.0 + Param.redshift);
  double dDdy  = Param.growth_dDdy(A);
  double dD2dy = Param.growth_dD2dy(A);
#endif

  lint np_saved = 0;
  Particle *p = particles->p;
  for(j = 0; j < particles->np_indomain; j++) {
    int ax;
    float x[3], v[3];
#ifdef PARTICLE_ID
    ulint id = (ulint) Param.P[j].ID;
#else
    ulint id = 0;
#endif
   
    if(UseCOLA){
#ifdef SCALEDEPENDENT
      // This assumes that dDdy acctually contains dD/dy instead DeltaD which we need for time-stepping
      // This is fine when we call it from the output-routine
      for(ax = 0; ax < 3; ax++){
        x[ax] = (float)(norm_pos*Param.P[j].Pos[ax]);
        v[ax] = (float)(norm_vel*(Param.P[j].Vel[ax] + (Param.P[j].dDdy[ax] + Param.P[j].dD2dy[ax])));
      }
#else
      for(ax = 0; ax < 3; ax++){
        x[ax] = (float)(norm_pos*Param.P[j].Pos[ax]);
        v[ax] = (float)(norm_vel*(Param.P[j].Vel[ax] + (Param.P[j].D[ax] * dDdy + Param.P[j].D2[ax] * dD2dy)));
      }
#endif
    } else {
      for(ax = 0; ax < 3; ax++){
        x[ax] = (float)(norm_pos*Param.P[j].Pos[ax]);
        v[ax] = (float)(norm_vel*Param.P[j].Vel[ax]);
      }
    }

    // We already know that all particles are in current slice
    x[0] -= edge_total_left;
    for(ax = 0; ax < 3; ax++){
      p[j].x[ax] = x[ax];
      p[j].v[ax] = v[ax];
    }
    p[j].id     = id;
    p[j].fof_id =  0;
    p[j].cll_id = -1;
  }
  np_saved = particles->np_indomain;
  mm_msg_printf(" Done translating [%li] particles\n", particles->np_indomain);

  // Check total number of particles
  unsigned long long np_saved_here = np_saved;
  unsigned long long np_global = 0;
  MPI_Reduce(&np_saved_here,&np_global,1,MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,0,MPI_COMM_WORLD);
  if(this_node == 0 && np_global != np_tot)
    mm_msg_abort(123,"Not all particles were saved %llu != %llu\n",np_global,np_tot);
  particles->np_indomain = np_saved;

  // Sort particles by x[0]
  mm_msg_printf(" Sorting particles and sharing buffer [%li] with neighbor nodes\n", particles->np_toleft);
  qsort(particles->p,particles->np_indomain,sizeof(Particle),compare_parts);

  // Communicate buffer
  MPI_Sendrecv(particles->p,particles->np_toleft,
      ParticleMPI,Param.i_node_left,tag,
      &(particles->p[particles->np_indomain]),particles->np_fromright,
      ParticleMPI,Param.i_node_right,tag,
      MPI_COMM_WORLD,&stat);

  // Add offset to received particles
  lint i_off = particles->np_indomain;
  for(i = 0; i < particles->np_fromright; i++){
    particles->p[i_off+i].x[0] += dx_domain;
  }

  mm_msg_printf(" Done\n\n");
}

