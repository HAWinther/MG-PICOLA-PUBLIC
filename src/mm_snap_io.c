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
#include <dirent.h>
#ifdef MATCHMAKER_USEFITS
#include <fitsio.h>
#endif
#include "mm_common.h"

/*
static float wrap_float(float x,float lb){
  if(x<0) return x+lb;
  else if(x>=lb) return x-lb;
  else return x;
}
*/

static FILE *my_fopen(char *fname,char *mode){
  FILE *f=fopen(fname,mode);
  if(f==NULL) {
    mm_msg_abort(123,"Couldn't open file %s\n",fname);
  }

  return f;
}

int compare_halo(const void *p1,const void *p2){
  const FoFHalo *f1=p1;
  const FoFHalo *f2=p2;

  if(f1->np>f2->np)
    return -1;
  else if(f1->np<f2->np)
    return 1;
  else
    return 0;
}

void write_halos_binary(char *fname,lint n_halos,long n_halos_total,FoFHalo *fhal){
  int bsize;
  FoFHeader head;
  FILE *fo=my_fopen(fname,"wb");

  head.n_halos_total=(long)n_halos_total;
  head.n_halos_here=(long)n_halos;
  if(Param.output_pernode==1)
    head.n_files=Param.n_nodes;
  else
    head.n_files=1;

  head.boxsize  = Param.boxsize;;
  head.redshift = Param.redshift;
  head.omega_m  = Param.omega_m;
  head.omega_l  = Param.omega_l;
  head.hubble   = Param.h;

  bsize=sizeof(FoFHeader);
  fwrite(&bsize,sizeof(int),1,fo);
  fwrite(&head,sizeof(FoFHeader),1,fo);
  fwrite(&bsize,sizeof(int),1,fo);

  bsize=n_halos*sizeof(FoFHalo);
  fwrite(&bsize,sizeof(int),1,fo);
  fwrite(fhal,sizeof(FoFHalo),n_halos,fo);
  fwrite(&bsize,sizeof(int),1,fo);
  fclose(fo);
}

void write_halos_ascii(char *fname,lint n_halos,FoFHalo *fhal)
{
  lint ii;
  FILE *fo;

  fo=my_fopen(fname,"w");
  if(Param.i_node==0)
    fprintf(fo,"#ID[1] NP[2] Mass[3] x_avg[4,5,6] x_rms[7,8,9] v_avg[10,11,12] v_rms[13,14,15] b[16] c[17] L[18,19,20] Ea[21,22,23] Eb[24,25,26] Ec[27,28,29] Units: Mass in [Msun/h] Pos are in [%E Mpc/h], Vel are comoving dx/dt in [km/s]\n", 1.0 / Param.norm_pos);
  for(ii=0;ii<n_halos;ii++) {
    fprintf(fo,"%ld %d %E ",(long)ii,fhal[ii].np,fhal[ii].m_halo);
    fprintf(fo,"%E %E %E ",
        fhal[ii].x_avg[0],fhal[ii].x_avg[1],fhal[ii].x_avg[2]);
    fprintf(fo,"%E %E %E ",
        fhal[ii].x_rms[0],fhal[ii].x_rms[1],fhal[ii].x_rms[2]);
    fprintf(fo,"%E %E %E ",
        fhal[ii].v_avg[0],fhal[ii].v_avg[1],fhal[ii].v_avg[2]);
    fprintf(fo,"%E %E %E ",
        fhal[ii].v_rms[0],fhal[ii].v_rms[1],fhal[ii].v_rms[2]);
    fprintf(fo,"%E %E ",
        fhal[ii].b,fhal[ii].c);
    fprintf(fo,"%E %E %E ",
        fhal[ii].lam[0],fhal[ii].lam[1],fhal[ii].lam[2]);
    fprintf(fo,"%E %E %E ",
        fhal[ii].ea[0],fhal[ii].ea[1],fhal[ii].ea[2]);
    fprintf(fo,"%E %E %E ",
        fhal[ii].eb[0],fhal[ii].eb[1],fhal[ii].eb[2]);
    fprintf(fo,"%E %E %E ",
        fhal[ii].ec[0],fhal[ii].ec[1],fhal[ii].ec[2]);
    fprintf(fo,"\n");
  }
  fclose(fo);
}

void write_halos_fits(char *fname,lint n_halos,FoFHalo *fhal){
#ifdef MATCHMAKER_USEFITS
  lint ii;
  fitsfile *fptr;
  lint *lic_write;
  int *ic_write;
  float *fc_write;
  int status=0;
  char *ttype[]={"ID","NP","MASS",
    "PX_CM","PY_CM","PZ_CM","PX_RMS","PY_RMS","PZ_RMS",
    "VX_CM","VY_CM","VZ_CM","VX_RMS","VY_RMS","VZ_RMS",
    "LX","LY","LZ","B","C",
    "EAX","EAY","EAZ","EBX","EBY","EBZ","ECX","ECY","ECZ"};
#ifdef MATCHMAKER_LONGINT
  char *tform[]={"1J","1K","1E",
    "1E","1E","1E","1E","1E","1E",
    "1E","1E","1E","1E","1E","1E",
    "1E","1E","1E","1E","1E",
    "1E","1E","1E","1E","1E","1E","1E","1E","1E"};
#else //MATCHMAKER_LONGINT
  char *tform[]={"1J","1J","1E",
    "1E","1E","1E","1E","1E","1E",
    "1E","1E","1E","1E","1E","1E",
    "1E","1E","1E","1E","1E",
    "1E","1E","1E","1E","1E","1E","1E","1E","1E"};
#endif //MATCHMAKER_LONGINT
  char *tunit[]={"NA","NA","M_SUN_H",
    "MPC_H","MPC_H","MPC_H","MPC_H","MPC_H","MPC_H",
    "KM_S","KM_S","KM_S","KM_S","KM_S","KM_S",
    "MPC_H_KM_S","MPC_H_KM_S","MPC_H_KM_S","NA","NA",
    "NA","NA","NA","NA","NA","NA","NA","NA","NA"};
  long n_halosl=(long)n_halos;


  fits_create_file(&fptr,fname,&status);
  fits_create_tbl(fptr,BINARY_TBL,0,29,ttype,tform,tunit,NULL,&status);

  lic_write = (lint *) my_malloc(n_halos*sizeof(lint));
  if(lic_write==NULL)
    mm_msg_abort(123,"Out of memory\n");
  for(ii=0;ii<n_halos;ii++)  //IDs
    lic_write[ii]=(lint)ii;
#ifdef MATCHMAKER_LONGINT
  fits_write_col(fptr,TLONG,1 ,1,1,n_halosl,lic_write,&status);
#else //MATCHMAKER_LONGINT
  fits_write_col(fptr,TINT,1 ,1,1,n_halosl,lic_write,&status);
#endif //MATCHMAKER_LONGINT
  my_free(lic_write);

  ic_write = (int *) my_malloc(n_halos*sizeof(int));
  if(ic_write==NULL)
    mm_msg_abort(123,"Out of memory\n");
  for(ii=0;ii<n_halos;ii++)  {
    ic_write[ii]=(int)(fhal[ii].np);
    if(fhal[ii].np==0)
      mm_msg_abort(123,"WTF\n");
  }
  fits_write_col(fptr,TINT,2 ,1,1,n_halosl,ic_write,&status);
  my_free(ic_write);

  fc_write = (float *) my_malloc(n_halos*sizeof(float));
  if(fc_write==NULL)
    mm_msg_abort(123,"Out of memory\n");
  for(ii=0;ii<n_halos;ii++) //Mass
    fc_write[ii]=(float)(fhal[ii].m_halo);
  fits_write_col(fptr,TFLOAT,3 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //x_CM
    fc_write[ii]=(float)(fhal[ii].x_avg[0]);
  fits_write_col(fptr,TFLOAT,4 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //y_CM
    fc_write[ii]=(float)(fhal[ii].x_avg[1]);
  fits_write_col(fptr,TFLOAT,5 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //z_CM
    fc_write[ii]=(float)(fhal[ii].x_avg[2]);
  fits_write_col(fptr,TFLOAT,6 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //x_RMS
    fc_write[ii]=(float)(fhal[ii].x_rms[0]);
  fits_write_col(fptr,TFLOAT,7 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //y_RMS
    fc_write[ii]=(float)(fhal[ii].x_rms[1]);
  fits_write_col(fptr,TFLOAT,8 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //z_RMS
    fc_write[ii]=(float)(fhal[ii].x_rms[2]);
  fits_write_col(fptr,TFLOAT,9 ,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //vx_CM
    fc_write[ii]=(float)(fhal[ii].v_avg[0]);
  fits_write_col(fptr,TFLOAT,10,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //vy_CM
    fc_write[ii]=(float)(fhal[ii].v_avg[1]);
  fits_write_col(fptr,TFLOAT,11,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //vz_CM
    fc_write[ii]=(float)(fhal[ii].v_avg[2]);
  fits_write_col(fptr,TFLOAT,12,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //vx_RMS
    fc_write[ii]=(float)(fhal[ii].v_rms[0]);
  fits_write_col(fptr,TFLOAT,13,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //vy_RMS
    fc_write[ii]=(float)(fhal[ii].v_rms[1]);
  fits_write_col(fptr,TFLOAT,14,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //vz_RMS
    fc_write[ii]=(float)(fhal[ii].v_rms[2]);
  fits_write_col(fptr,TFLOAT,15,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //Lx
    fc_write[ii]=(float)(fhal[ii].lam[0]);
  fits_write_col(fptr,TFLOAT,16,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //Ly
    fc_write[ii]=(float)(fhal[ii].lam[1]);
  fits_write_col(fptr,TFLOAT,17,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //Lz
    fc_write[ii]=(float)(fhal[ii].lam[2]);
  fits_write_col(fptr,TFLOAT,18,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //B
    fc_write[ii]=(float)(fhal[ii].b);
  fits_write_col(fptr,TFLOAT,19,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //C
    fc_write[ii]=(float)(fhal[ii].c);
  fits_write_col(fptr,TFLOAT,20,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //EAx
    fc_write[ii]=(float)(fhal[ii].ea[0]);
  fits_write_col(fptr,TFLOAT,21,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //EAy
    fc_write[ii]=(float)(fhal[ii].ea[1]);
  fits_write_col(fptr,TFLOAT,22,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //EAz
    fc_write[ii]=(float)(fhal[ii].ea[2]);
  fits_write_col(fptr,TFLOAT,23,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //EBx
    fc_write[ii]=(float)(fhal[ii].eb[0]);
  fits_write_col(fptr,TFLOAT,24,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //EBy
    fc_write[ii]=(float)(fhal[ii].eb[1]);
  fits_write_col(fptr,TFLOAT,25,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //EBz
    fc_write[ii]=(float)(fhal[ii].eb[2]);
  fits_write_col(fptr,TFLOAT,26,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //ECx
    fc_write[ii]=(float)(fhal[ii].ec[0]);
  fits_write_col(fptr,TFLOAT,27,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //ECy
    fc_write[ii]=(float)(fhal[ii].ec[1]);
  fits_write_col(fptr,TFLOAT,28,1,1,n_halosl,fc_write,&status);
  for(ii=0;ii<n_halos;ii++) //ECz
    fc_write[ii]=(float)(fhal[ii].ec[2]);
  fits_write_col(fptr,TFLOAT,29,1,1,n_halosl,fc_write,&status);
  my_free(fc_write);
  fits_close_file(fptr,&status);
#else
  mm_msg_abort(123, "Error: MatchMaker not compiled with FITS! Try a different output-format!\n");
#endif
}

void write_halos_multi(char *prefix,lint n_halos,FoFHalo *fhal,long n_halos_total){
  char fname[256];

  if(Param.output_format==0) {
    sprintf(fname,"%s.txt",prefix);
    mm_msg_printf(" Writing to file %s\n",fname);
    write_halos_ascii(fname,n_halos,fhal);
  } else if(Param.output_format==1) {
    sprintf(fname,"%s.fits",prefix);
    mm_msg_printf(" Writing to file %s\n",fname);
    write_halos_fits(fname,n_halos,fhal);
  } else if(Param.output_format==2) {
    sprintf(fname,"%s.dat",prefix);
    mm_msg_printf(" Writing to file %s\n",fname);
    write_halos_binary(fname,n_halos,n_halos_total,fhal);
  }
}

static void write_halos_pernode(lint n_halos,FoFHalo *fhal){
  char prefix[256];
  long n_halos_total;
  long n_halos_here=n_halos;

  MPI_Allreduce(&n_halos_here,&n_halos_total,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  mm_msg_printf("%ld halos found in total\n",n_halos_total);

  //Sort by mass
  qsort(fhal,n_halos,sizeof(FoFHalo),compare_halo);

  sprintf(prefix,"%s.%d",Param.output_prefix,Param.i_node);
  write_halos_multi(prefix,n_halos,fhal,n_halos_total);

  my_free(fhal);
}

static void write_halos_root(lint n_halos,FoFHalo *fhal){
  //Gather halos
  int ii;
  int n_halos_total=0;
  int n_halos_here=(int)n_halos;
  int *n_halos_array=NULL;
  int *n_disp_array=NULL;
  FoFHalo *fh_tot=NULL;

  if(Param.i_node==0) {
    n_halos_array = my_malloc(Param.n_nodes*sizeof(int));
    if(n_halos_array==NULL)
      mm_msg_abort(123,"Out of memory\n");
    n_disp_array = my_malloc(Param.n_nodes*sizeof(int));
    if(n_disp_array==NULL)
      mm_msg_abort(123,"Out of memory\n");
  }

  MPI_Gather(&n_halos_here,1,MPI_INT,
      n_halos_array,1,MPI_INT,
      0,MPI_COMM_WORLD);

  if(Param.i_node==0) {
    for(ii=0;ii<Param.n_nodes;ii++)
      n_halos_total+=n_halos_array[ii];
    if(n_halos_total<=0){
      mm_msg_abort(123,"No haloes found\n");
    } else {
      mm_msg_printf(" %ld halos found in total with %d particles of more\n",
          (long)n_halos_total,Param.np_min);
    }
    //Calculate array of displacements
    n_disp_array[0]=0;
    for(ii=1;ii<Param.n_nodes;ii++)
      n_disp_array[ii]=n_disp_array[ii-1]+n_halos_array[ii-1];

    fh_tot = my_malloc(n_halos_total*sizeof(FoFHalo));
    if(fh_tot==NULL)
      mm_msg_abort(123,"Failed to allocate memory for all halos\n");
  }
  MPI_Gatherv(fhal,n_halos_here,HaloMPI,
      fh_tot,n_halos_array,n_disp_array,HaloMPI,
      0,MPI_COMM_WORLD);
  my_free(fhal);

  if(Param.i_node==0) {
    char prefix[256];

    my_free(n_halos_array);
    my_free(n_disp_array);

    //Sort by mass
    qsort(fh_tot,n_halos_total,sizeof(FoFHalo),compare_halo);

    sprintf(prefix,"%s",Param.output_prefix);
    write_halos_multi(prefix,n_halos_total,fh_tot,(long)n_halos_total);

    my_free(fh_tot);
  }

  mm_msg_printf(" Done\n\n");
}

void write_halos(lint n_halos,FoFHalo *fhal){
  if(Param.output_pernode || n_halos == 0){
    write_halos_pernode(n_halos,fhal);
  } else {
    write_halos_root(n_halos,fhal);
  }
}
