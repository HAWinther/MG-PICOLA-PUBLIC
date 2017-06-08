#ifndef MATCHMAKER_COMMON_H
#define MATCHMAKER_COMMON_H
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

#include <mpi.h>
#include "wrappers.h"

#define MATCHMAKER_LONGINT

#define MATCHMAKER_LONGIDS

#define MATCHMAKER_MASS_FACTOR 1.0E10

#define N_IDS_NEW_MAX 1000000

//#define MATCHMAKER_DEBUG

//#define MATCHMALER_CORRECT_NP

#ifdef MATCHMAKER_LONGINT
typedef long lint;
typedef unsigned long long ulint;
#else
typedef int lint;
typedef unsigned int ulint;
#endif

typedef struct {
  char output_prefix[256];   // Output prefix
  int output_pernode;
  double dx_extra;           // Buffer for parallelization (in Mpc/h)
  double b_fof;              // Percolation fraction
  int np_min;                // Minimum number of particles
  int output_format;         // Output format
  
  char OutputDir[500];       
  char FileBase[500];

  unsigned long long n_part; // Total number of particles
  int n_part_1d;             // Total number of particles in 1D
  int NumPart;               // Number of particles in domain
  double boxsize;            // Box size
  double omega_m;            // Omega_M
  double omega_l;            // Omega_Lambda
  double redshift;           // Redshift
  double norm_vel;           // From PICOLA v to peculiar v in km/s
  double norm_pos;           // From PICOLA pos to pos in Mpc/h
  double h;                  // Hubble parameter
  double mp;                 // Particle mass
  
  const struct part_data *P; // PICOLA particles
  
  int Local_p_start;   // Index for where the PICOLA particle-grid starts

  // Growth-functions of a
  double (*growth_dDdy)(double);
  double (*growth_dD2dy)(double);

  // MPI parameters
  int n_nodes;
  int i_node;
  int i_node_left;
  int i_node_right;
  float x_offset;
  float x_offset_right;
  float dx_domain;
} Parameters;
Parameters Param;

typedef struct {
  float x[3];
  float v[3];
#ifdef MATCHMAKER_LONGIDS
  unsigned long long id;
#else  // MATCHMAKER_LONGIDS
  unsigned int id;
#endif // MATCHMAKER_LONGIDS
  lint fof_id;
  lint cll_id;
} Particle;

typedef struct {
  Particle* p;
  Particle* p_back;
  lint np_allocated; // Total number of particles allocated
  lint np_local;     // Total number of particles saved
  lint np_indomain;  // Total number of particles found in domain
  lint np_toleft;    // Number of particles that will be passed on to the node on the left
  lint np_centre;    // Number of particles that won't be passed to nor from any other node
  lint np_fromright; // Number of particles passed from the right node
  unsigned long long np_total;
  float np_average;
} Particles;

typedef struct {
  int np;
  float m_halo;
  float x_avg[3];
  float x_rms[3];
  float v_avg[3];
  float v_rms[3];
  float lam[3];
  float b;
  float c;
  float ea[3];
  float eb[3];
  float ec[3];
} FoFHalo;

typedef struct {
  long n_halos_total;
  long n_halos_here;
  int n_files;
  float boxsize;
  float redshift;
  float omega_m;
  float omega_l;
  float hubble;
  char fill[256-40];
} FoFHeader;

MPI_Datatype ParticleMPI;
MPI_Datatype HaloMPI;

FoFHalo *fof_get_halos(lint *n_halos_out,Particles *particles);
void picola_to_matchmaker_particles(Particles *particles);
void write_halos(lint n_halos,FoFHalo *fhal);
void mm_msg_init();
void mm_msg_printf(const char *fmt, ...);
void mm_msg_abort(const int errret,const char *fmt, ...);

#endif
