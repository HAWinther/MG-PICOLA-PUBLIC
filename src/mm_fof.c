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
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <dirent.h>
#include "mm_common.h"

static lint Ngrid[3];
static lint Ngrid_tot;
static float Ac[3];
static float I_ac[3];
static float Dfof;
static float D2fof;
static float Lbox;
static float Lbox_half;

static lint *ids_new,*ids_old,*ids_single;

typedef struct {
  double eval;
  double evec[3];
} Esys;

typedef struct {
  lint np_in;
  lint *ids_in;
  lint np_out;
  lint *ids_out;
} Cell;

typedef struct {
  int np;
  lint *ids;
} FoFGroup;

static int compare_evals(const void *p1,const void *p2){
  Esys *e1=(Esys *)p1;
  Esys *e2=(Esys *)p2;

  return (e1->eval<e2->eval) - (e1->eval>e2->eval);
}

static int compare_fof(const void *p1,const void *p2){
  const FoFGroup *f1=p1;
  const FoFGroup *f2=p2;

  if(f1->np>f2->np)
    return -1;
  else if(f1->np<f2->np)
    return 1;
  else
    return 0;
}

void init_fof(void){
  int ax;
  float ipd=Param.boxsize/pow((double)(Param.n_part),1./3.);
  float dx[3];

  dx[0]=Param.dx_domain+Param.dx_extra;
  dx[1]=Param.boxsize;
  dx[2]=Param.boxsize;

  for(ax=0;ax<3;ax++) {
    Ngrid[ax]=(lint)(dx[ax]/ipd)+1;
    Ac[ax]=dx[ax]/Ngrid[ax];
    I_ac[ax]=1./Ac[ax];
  }

  Ngrid_tot=Ngrid[0]*Ngrid[1]*Ngrid[2];
  Dfof=ipd*Param.b_fof;
  D2fof=Dfof*Dfof;
  Lbox=(float)(Param.boxsize);
  Lbox_half=(float)(Param.boxsize/2);
#ifdef MATCHMAKER_DEBUG
  mm_msg_printf(" Grid size = %d (%d x %d x %d)\n",
      Ngrid_tot,Ngrid[0],Ngrid[1],Ngrid[2]);
  mm_msg_printf(" Cell sizes : %lf x %lf x %lf\n",
      Ac[0],Ac[1],Ac[2]);
  mm_msg_printf(" FoF separation : %lf\n",Dfof);
#endif //MATCHMAKER_DEBUG
}

static int get_neighbors(lint ip0,Cell *cll,Particle *p,lint fof_id)
{
  lint ii;
  Particle *p0=&(p[ip0]);
  Cell c0=cll[p0->cll_id];
  float *x0=p0->x;
  int n_new_neighbors=0;

  //Couple with particles in the same cell
  for(ii=0;ii<c0.np_in;ii++) {
    lint id=c0.ids_in[ii];

    if(id!=ip0) {
      if(p[id].fof_id==0) { //DEBUG THIS
        float *x1=p[id].x;
        float dx;

        dx=fabs(x0[0]-x1[0]);
        if(dx<=Dfof) {
          float dy,dz,d2;
          dy=fabs(x0[1]-x1[1]);
          dz=fabs(x0[2]-x1[2]);
          if(dy>Lbox_half) dy=Lbox-dy;
          if(dz>Lbox_half) dz=Lbox-dz;
          d2=dx*dx+dy*dy+dz*dz;
          if(d2<=D2fof) {
            if(p0->fof_id==0) //New halo!
              p0->fof_id=fof_id;
            p[id].fof_id=fof_id;
            ids_single[n_new_neighbors]=id;
            n_new_neighbors++;
            if(n_new_neighbors>=N_IDS_NEW_MAX)
              mm_msg_abort(123,"Enlarge N_IDS_NEW_MAX\n");
          }
        }
      }
    }
  }

  //Couple with particles in the same cell
  for(ii=0;ii<c0.np_out;ii++) {
    lint id=c0.ids_out[ii];

    if(id!=ip0) {
      if(p[id].fof_id==0) { //DEBUG THIS
        float dx;
        float *x1=p[id].x;

        dx=fabs(x0[0]-x1[0]);
        if(dx<=Dfof) {
          float dy,dz,d2;
          dy=fabs(x0[1]-x1[1]);
          dz=fabs(x0[2]-x1[2]);
          if(dy>Lbox_half) dy=Lbox-dy;
          if(dz>Lbox_half) dz=Lbox-dz;
          d2=dx*dx+dy*dy+dz*dz;
          if(d2<=D2fof) {
            if(p0->fof_id==0) //New halo!
              p0->fof_id=fof_id;
            p[id].fof_id=fof_id;
            ids_single[n_new_neighbors]=id;
            n_new_neighbors++;
            if(n_new_neighbors>=N_IDS_NEW_MAX)
              mm_msg_abort(123,"Enlarge N_IDS_NEW_MAX\n");
          }
        }
      }
    }
    else
      mm_msg_abort(123,"This shouldn't have happened\n");
  }

  return n_new_neighbors;
}

static Cell *gather_particles_in_cells(Particles *particles) {
  lint i;

  if(Param.i_node == 0){
    printf(" Allocated memory per task: [%6.3f] MB (particles) [%6.3f] MB (grid)\n", (particles->np_toleft+particles->np_local)*sizeof(Particle)/1e6, Ngrid_tot*sizeof(Cell)/1e6);
  }

  Cell *cll = my_malloc(Ngrid_tot*sizeof(Cell));
  if(cll==NULL)
    mm_msg_abort(123,"Failed to allocate memory for NGBS cells\n");

  for(i=0;i<Ngrid_tot;i++) {
    cll[i].np_in=0;
    cll[i].np_out=0;
  }

  Particle *p=particles->p;
  for(i=0;i<particles->np_local;i++) {
    int ax;
    lint ic[3],ic_close[3];
    lint icell;

    for(ax=0;ax<3;ax++) {
      float x=p[i].x[ax];
      float ix=x*I_ac[ax];
      lint ic_in=(lint)ix;
      lint ic_cl=(lint)(ix+0.5);
      float dist_close=fabs(ic_cl*Ac[ax]-x);

      if(ic_in>=Ngrid[ax]) {
        if(ax==0)
          ic[ax]=Ngrid[ax]-1;
        else
          ic[ax]=ic_in-Ngrid[ax];
      }
      else if(ic_in<0) { 
        if(ax==0)
          ic[ax]=0;
        else
          ic[ax]=ic_in+Ngrid[ax];
      }
      else ic[ax]=ic_in;

      if(dist_close<=Dfof) {
        ic_close[ax]=ic[ax]-1+2*(ic_cl-ic_in);
        if(ic_close[ax]>=Ngrid[ax]) {
          if(ax==0)
            ic_close[ax]=-1;
          else
            ic_close[ax]-=Ngrid[ax];
        }
        else if(ic_close[ax]<0) {
          if(ax==0)
            ic_close[ax]=-1;
          else
            ic_close[ax]+=Ngrid[ax];
        }
      }
      else
        ic_close[ax]=-1;
    }
    icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
    cll[icell].np_in++;
    p[i].cll_id=icell;
    p[i].fof_id=0;

    if(ic_close[0]>=0) {
      icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
      cll[icell].np_out++;
      if(ic_close[1]>=0) {
        icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
        cll[icell].np_out++;
        icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
        cll[icell].np_out++;
        if(ic_close[2]>=0) {
          icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
          cll[icell].np_out++;
          icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
          cll[icell].np_out++;
          icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
          cll[icell].np_out++;
          icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
          cll[icell].np_out++;
        }
      }
      else if(ic_close[2]>=0) {
        icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
        cll[icell].np_out++;
        icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
        cll[icell].np_out++;
      }
    }
    else if(ic_close[1]>=0) {
      icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
      cll[icell].np_out++;
      if(ic_close[2]>=0) {
        icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
        cll[icell].np_out++;
        icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
        cll[icell].np_out++;
      }
    }
    else if(ic_close[2]>=0) {
      icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
      cll[icell].np_out++;
    }
  }

  // NB: We don't use [my_malloc] here as this is too slow
  for(i=0;i<Ngrid_tot;i++) {
    lint np_in=cll[i].np_in;
    lint np_out=cll[i].np_out;
    if(np_in>0) {
      cll[i].ids_in=malloc(np_in*sizeof(lint));
      if(cll[i].ids_in==NULL)
        mm_msg_abort(123,"Failed to allocate ids in cells\n");
    }
    if(np_out>0) {
      cll[i].ids_out=malloc(np_out*sizeof(lint));
      if(cll[i].ids_out==NULL)
        mm_msg_abort(123,"Failed to allocate ids in cells\n");
    }
    cll[i].np_in=0;
    cll[i].np_out=0;
  }

  for(i=0;i<particles->np_local;i++) {
    int ax;
    lint ic[3],ic_close[3];
    lint icell;

    for(ax=0;ax<3;ax++) {
      float x=p[i].x[ax];
      float ix=x*I_ac[ax];
      lint ic_in=(lint)ix;
      lint ic_cl=(lint)(ix+0.5);
      float dist_close=fabs(ic_cl*Ac[ax]-x);

      if(ic_in>=Ngrid[ax]) {
        if(ax==0)
          ic[ax]=Ngrid[ax]-1;
        else
          ic[ax]=ic_in-Ngrid[ax];
      }
      else if(ic_in<0) { 
        if(ax==0)
          ic[ax]=0;
        else
          ic[ax]=ic_in+Ngrid[ax];
      }
      else ic[ax]=ic_in;

      if(dist_close<=Dfof) {
        ic_close[ax]=ic[ax]-1+2*(ic_cl-ic_in);
        if(ic_close[ax]>=Ngrid[ax]) {
          if(ax==0)
            ic_close[ax]=-1;
          else
            ic_close[ax]-=Ngrid[ax];
        }
        else if(ic_close[ax]<0) {
          if(ax==0)
            ic_close[ax]=-1;
          else
            ic_close[ax]+=Ngrid[ax];
        }
      }
      else
        ic_close[ax]=-1;
    }
    icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
    if(icell!=p[i].cll_id)
      mm_msg_abort(123,"Error in cell assignment, node %d\n",Param.i_node);
    cll[icell].ids_in[cll[icell].np_in++]=i;

    if(ic_close[0]>=0) {
      icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
      cll[icell].ids_out[cll[icell].np_out++]=i;
      if(ic_close[1]>=0) {
        icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
        cll[icell].ids_out[cll[icell].np_out++]=i;
        icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
        cll[icell].ids_out[cll[icell].np_out++]=i;
        if(ic_close[2]>=0) {
          icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
          cll[icell].ids_out[cll[icell].np_out++]=i;
          icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
          cll[icell].ids_out[cll[icell].np_out++]=i;
          icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
          cll[icell].ids_out[cll[icell].np_out++]=i;
          icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
          cll[icell].ids_out[cll[icell].np_out++]=i;
        }
      }
      else if(ic_close[2]>=0) {
        icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
        cll[icell].ids_out[cll[icell].np_out++]=i;
        icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
        cll[icell].ids_out[cll[icell].np_out++]=i;
      }
    }
    else if(ic_close[1]>=0) {
      icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
      cll[icell].ids_out[cll[icell].np_out++]=i;
      if(ic_close[2]>=0) {
        icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
        cll[icell].ids_out[cll[icell].np_out++]=i;
        icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
        cll[icell].ids_out[cll[icell].np_out++]=i;
      }
    }
    else if(ic_close[2]>=0) {
      icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
      cll[icell].ids_out[cll[icell].np_out++]=i;
    }
  }

  return cll;
}

static FoFGroup *assign_particles_to_fof(Particles *particles,Cell *cll,lint *n_fof_out){
  lint i;
  lint n_fof=1;
  Particle *p=particles->p;

  ids_new=(lint *)my_malloc(N_IDS_NEW_MAX*sizeof(lint));
  if(ids_new==NULL)
    mm_msg_abort(123,"Out of memory\n");
  ids_old=(lint *)my_malloc(N_IDS_NEW_MAX*sizeof(lint));
  if(ids_old==NULL)
    mm_msg_abort(123,"Out of memory\n");
  ids_single=(lint *)my_malloc(N_IDS_NEW_MAX*sizeof(lint));
  if(ids_single==NULL)
    mm_msg_abort(123,"Out of memory\n");
  //Do FoF
  for(i=0;i<particles->np_indomain;i++) {
    if(p[i].fof_id==0) {
      int n_old=get_neighbors(i,cll,p,n_fof);
      if(n_old>0) {
        memcpy(ids_old,ids_single,n_old*sizeof(lint));
        while(n_old>0) {
          int ii;
          int n_new=0;
          for(ii=0;ii<n_old;ii++) {
            lint id=ids_old[ii];
            int n_this=get_neighbors(id,cll,p,n_fof);
            memcpy(&(ids_new[n_new]),ids_single,n_this*sizeof(lint));
            n_new+=n_this;
          }
          n_old=n_new;
          memcpy(ids_old,ids_new,n_new*sizeof(lint));
        }
      }
      if(p[i].fof_id==n_fof) //Halo found
        n_fof++;
    }
  }
  n_fof--;
  my_free(ids_single);
  my_free(ids_new);
  my_free(ids_old);

  // Free memory from cells
  for(i=0;i<Ngrid_tot;i++) {
    if(cll[i].np_in>0)
      free(cll[i].ids_in);
    if(cll[i].np_out>0)
      free(cll[i].ids_out);
  }
  my_free(cll);

  //Allocate FoF groups
  FoFGroup *fof=my_malloc(n_fof*sizeof(FoFGroup));
  if(fof==NULL)
    mm_msg_abort(123,"Unable to allocate memory for FoF groups\n");
  for(i=0;i<n_fof;i++)
    fof[i].np=0;

  //Compare particles
  //Receive buffer particles from left and send to right
  int tag=200;
  MPI_Status stat;
  MPI_Sendrecv(&(particles->p[particles->np_indomain]),particles->np_fromright,
      ParticleMPI,Param.i_node_right,tag,
      particles->p_back,particles->np_toleft,
      ParticleMPI,Param.i_node_left,tag,
      MPI_COMM_WORLD,&stat);

  //Resolve redundancies
  Particle *p_fromleft=particles->p_back;
  for(i=0;i<particles->np_toleft;i++) {
    if(p_fromleft[i].fof_id!=0)
      p[i].fof_id=0;
  }

  //Allocate id arrays in halos
  for(i=0;i<particles->np_local;i++) {
    if(p[i].fof_id!=0)
      fof[p[i].fof_id-1].np++;
  }

  lint n_fof_true=0;
  for(i=0;i<n_fof;i++) {
    int np=fof[i].np;
    if(np>0) {
      n_fof_true++;
      fof[i].ids=malloc(np*sizeof(lint));
      if(fof[i].ids==NULL)
        mm_msg_abort(123,"Failed to allocate memory for particle ids in FoF groups\n");
      fof[i].np=0;
    }
  }
#ifdef MATCHMAKER_DEBUG
  printf("In node %d there were initially %ld FoF groups "
      "but only %ld after talking to node %d. Solved %ld redundancies\n",
      Param.i_node,(long)n_fof,(long)n_fof_true,
      Param.i_node_left,(long)(n_fof-n_fof_true));
#endif //MATCHMAKER_DEBUG  

  //Assign particle ids to FoF groups
  for(i=0;i<particles->np_local;i++) {
    if(p[i].fof_id!=0) {
      lint id=i;
      lint i_fof=p[i].fof_id-1;
      fof[i_fof].ids[fof[i_fof].np]=id;
      fof[i_fof].np++;
    }
  }

  //Sort by number of particles
  qsort(fof,n_fof,sizeof(FoFGroup),compare_fof);

  *n_fof_out=n_fof_true;

  return fof;
}

static FoFHalo *get_halos(lint *n_halos_out,lint n_fof,FoFGroup *fg,Particles *particles)
{
  lint i,n_halos;
  Particle *p=particles->p;
  FoFHalo *fh;
  gsl_matrix *inertia=gsl_matrix_alloc(3,3);
  gsl_matrix *evec=gsl_matrix_alloc(3,3);
  gsl_vector *eval=gsl_vector_alloc(3);
  gsl_eigen_symmv_workspace *ws=gsl_eigen_symmv_alloc(3);
  Esys leig[3];

  lint np_current=fg[0].np;
  n_halos=0;
  for(i=0;i<n_fof;i++) {
    lint np=fg[i].np;
    if(np>np_current)
      mm_msg_abort(123,"Ordering seems to be wrong!\n");
    if(np>=Param.np_min)
      n_halos++;
    np_current=np;
  }
  if(n_halos==0) {
#ifdef MATCHMAKER_DEBUG
    printf("Node %d couldn't find any halos above np_min=%d\n",Param.i_node,Param.np_min);
#endif //MATCHMAKER_DEBUG  
    *n_halos_out=0;
    return NULL;
  }

#ifdef MATCHMAKER_DEBUG
  printf("In node %d there were initially %ld FoF groups "
      "but only %ld of them have >%d particles\n",
      Param.i_node,(long)n_fof,(long)n_halos,Param.np_min);
#endif //MATCHMAKER_DEBUG  

  fh=my_malloc(n_halos*sizeof(FoFHalo));
  if(fh==NULL)
    mm_msg_abort(123,"Unable to allocate memory for halos\n");  

  double mass_particle=MATCHMAKER_MASS_FACTOR*Param.mp;
  for(i=0;i<n_halos;i++) {
    lint j;
    double inertia_here[9];
    lint np=fg[i].np;
    double x_avg[3],v_avg[3],x_rms[3],v_rms[3],lam[3];

    if(np<2)
      printf("WTF\n");

    fh[i].np=np;
#ifdef MATCHMAKER_CORRECT_NP
    fh[i].m_halo=np*(1-pow((double)np,-0.6))*mass_particle;
#else //MATCHMAKER_CORRECT_NP
    fh[i].m_halo=np*mass_particle;
#endif //MATCHMAKER_CORRECT_NP
    for(j=0;j<3;j++) {
      int k;
      x_avg[j]=0;
      v_avg[j]=0;
      x_rms[j]=0;
      v_rms[j]=0;
      lam[j]=0;
      fh[i].x_avg[j]=0;
      fh[i].x_rms[j]=0;
      fh[i].v_avg[j]=0;
      fh[i].v_rms[j]=0;
      fh[i].lam[j]=0;
      for(k=0;k<3;k++)
        inertia_here[k+3*j]=0.;
    }

    //Compute center of mass position and velocity
    for(j=0;j<np;j++) {
      int ax;
      lint ip=fg[i].ids[j];
      for(ax=0;ax<3;ax++) {
        double x=p[ip].x[ax];
        double v=p[ip].v[ax];

        if(j>0) {
          //If the distance to the current center of mass is larger
          //than L/2, wrap around
          double cm_here=x_avg[ax]/j;
          if(2*fabs(x-cm_here)>Param.boxsize) {
            if(2*x>Param.boxsize)
              x-=Param.boxsize;
            else 
              x+=Param.boxsize;
          }
        }
        x_avg[ax]+=x;
        v_avg[ax]+=v;
      }
    }
    for(j=0;j<3;j++) {
      fh[i].x_avg[j]=x_avg[j]/np;
      fh[i].v_avg[j]=v_avg[j]/np;
    }

    //Compute relative quantities
    for(j=0;j<np;j++) {
      int ax;
      double dx[3];
      double dv[3];
      lint ip=fg[i].ids[j];

      for(ax=0;ax<3;ax++) {
        double x=p[ip].x[ax];
        double v=p[ip].v[ax];

        if(2*fabs(x-fh[i].x_avg[ax])>Param.boxsize) {
          //Wrap around if the distance to CM is larger than L/2
          if(2*x>Param.boxsize)
            x-=Param.boxsize;
          else
            x+=Param.boxsize;
        }
        dx[ax]=x-fh[i].x_avg[ax];
        dv[ax]=v-fh[i].v_avg[ax];
      }

      //RMS
      for(ax=0;ax<3;ax++) {
        x_rms[ax]+=dx[ax]*dx[ax];
        v_rms[ax]+=dv[ax]*dv[ax];
      }

      //Inertia tensor;
      for(ax=0;ax<3;ax++) {
        int ax2;
        for(ax2=0;ax2<3;ax2++)
          inertia_here[ax2+3*ax]+=dx[ax]*dx[ax2];
      }

      //Angular momentum
      lam[0]+=dx[1]*dv[2]-dx[2]*dv[1];
      lam[1]+=dx[2]*dv[0]-dx[0]*dv[2];
      lam[2]+=dx[0]*dv[1]-dx[1]*dv[0];
    }


    //Diagonalize inertia tensor
    for(j=0;j<3;j++) {
      int k;
      for(k=0;k<3;k++)
        gsl_matrix_set(inertia,j,k,inertia_here[k+3*j]);
    }
    gsl_eigen_symmv(inertia,eval,evec,ws);
    for(j=0;j<3;j++) {
      leig[j].eval=gsl_vector_get(eval,j);
      leig[j].evec[0]=gsl_matrix_get(evec,0,j);
      leig[j].evec[1]=gsl_matrix_get(evec,1,j);
      leig[j].evec[2]=gsl_matrix_get(evec,2,j);
    }
    qsort(leig,3,sizeof(Esys),compare_evals);

    //Add evals and evecs
    for(j=0;j<3;j++) {
      int k;
      if(leig[0].eval<=0) {
        fh[i].b=0;
        fh[i].c=0;
        for(k=0;k<3;k++) {
          fh[i].ea[k]=0;
          fh[i].eb[k]=0;
          fh[i].ec[k]=0;
        }
      }
      else {
        fh[i].b=leig[1].eval/leig[0].eval;
        fh[i].c=leig[2].eval/leig[0].eval;
        for(k=0;k<3;k++) {
          fh[i].ea[k]=leig[0].evec[k];
          fh[i].eb[k]=leig[1].evec[k];
          fh[i].ec[k]=leig[2].evec[k];
        }
      }
    }

    for(j=0;j<3;j++) {
      //Normalize
      fh[i].x_rms[j]=sqrt(x_rms[j]/np);
      fh[i].v_rms[j]=sqrt(v_rms[j]/np);
      fh[i].lam[j]=lam[j];

      //Wrap CM
      if(fh[i].x_avg[j]<0) //Wrap CM
        fh[i].x_avg[j]+=Param.boxsize;
      else if(fh[i].x_avg[j]>=Param.boxsize)
        fh[i].x_avg[j]-=Param.boxsize;
    }
    fh[i].x_avg[0]+=Param.x_offset;
  }

  for(i=0;i<n_fof;i++) {
    if(fg[i].np>0)
      free(fg[i].ids);
  }
  my_free(fg);

  gsl_matrix_free(inertia);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(ws);

  *n_halos_out = n_halos;
  return fh;
}

FoFHalo *fof_get_halos(lint *n_halos_out,Particles *particles){
  FoFGroup *fof_g;
  FoFHalo *fof_h;
  Cell *cll;
  lint n_fof,n_halos;

  init_fof();
  mm_msg_printf(" Gathering particles in cells for NB searching\n");
  cll = gather_particles_in_cells(particles);
  mm_msg_printf(" Matchmaking\n");
  fof_g = assign_particles_to_fof(particles,cll,&n_fof);
  mm_msg_printf(" Computing halo properties\n");
  fof_h = get_halos(&n_halos,n_fof,fof_g,particles);
  mm_msg_printf(" Done\n\n");
  *n_halos_out = n_halos;

  return fof_h;
}
