//==========================================================================//
//  Copyright (c) 2013       Cullan Howlett & Marc Manera,                  //
//                           Institute of Cosmology and Gravitation,        //
//                           University of Portsmouth.                      //
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

//==========================================================================//
// This file contains routines to read in and check the run parameters file //
//==========================================================================//

#include "vars.h"
#include "proto.h"

//============================================================================
// Read in the list of output redshifts and the number of steps between outputs
//============================================================================
void read_outputs(void) {

  FILE * fd;
  char buf[500];

  Noutputs = 0;
  if((fd = fopen(OutputRedshiftFile, "r"))) {
    fflush(stdout);
    while(fgets(buf,500,fd)) Noutputs++;
    fclose(fd);

    OutputList = (struct Outputs *) my_malloc(Noutputs*sizeof(struct Outputs));

    Noutputs = 0;
    fd = fopen(OutputRedshiftFile, "r");
    fflush(stdout);
    while(fgets(buf,500,fd)) {
      int nsteps;
      double red;
      if(buf[0] == '%') continue;
      if(sscanf(buf, "%lf, %d", &red, &nsteps) > 0) {
        if(sscanf(buf, "%lf, %d", &red, &nsteps) != 2) {
          if(ThisTask == 0) fprintf(stdout,"\nERROR: Line in Output Redshift File '%s' is in incorrect format.\n\n",buf);
          fclose(fd);
          FatalError((char *)"read_param.c Output Redshift File has wrong format");
        }
        if (nsteps <= 0) {
          if ((Init_Redshift-red)/Init_Redshift > 1.0E-6) {
            if(ThisTask == 0) fprintf(stdout,"\nERROR: I read a value for nsteps of <= 0 up to redshift %lf.\n\n", red);
            fclose(fd);
            FatalError((char *)"read_param.c nsteps <= 0");
          }
        }
        OutputList[Noutputs].Nsteps = nsteps;
        OutputList[Noutputs].Redshift = red;
        Noutputs++;
      }
    }
    fclose(fd);
  } else {
    if(ThisTask == 0) fprintf(stdout,"\nERROR: Output Redshift File '%s' not found.\n\n",OutputRedshiftFile);
    FatalError((char *)"read_param.c Output Redshift File not found");
  }

  if (Noutputs == 0) {
    if(ThisTask == 0) fprintf(stdout,"\nERROR: Found no output redshifts in file '%s'. Surely this is accidental?.\n\n",OutputRedshiftFile);
    FatalError((char *)"read_param.c Noutputs == 0");
  }

  // Sort the output list via the redshifts, in descending order (in case it already isn't)
  qsort(OutputList,Noutputs,sizeof(struct Outputs),sort_redshift);

  if (OutputList[0].Redshift > Init_Redshift) {
    if(ThisTask == 0) fprintf(stdout,"\nERROR: The highest output redshift (%lf) is greater than the initial redshift (%lf).\nzn", OutputList[0].Redshift, Init_Redshift);
    FatalError((char *)"read_param.c OutputList[0].Redshift > Init_Redshift");
  }

}

//=========================================================================
// The comparison function to sort the output redshifts in descending order
//=========================================================================
int sort_redshift(const void * Item1, const void * Item2) {
  struct Outputs * Output1 = (struct Outputs *)Item1;
  struct Outputs * Output2 = (struct Outputs *)Item2;
  if((*Output1).Redshift > (*Output2).Redshift) return -1;
  if((*Output1).Redshift < (*Output2).Redshift) return 1;
  if(ThisTask == 0) fprintf(stdout,"\nERROR: Duplicate output redshift (%lf) in Output Redshift File.\n\n", (*Output1).Redshift);
  FatalError((char *)"read_param.c sort_redshift, Duplicate output redshift");
  return 0;
}

//============================================================================
// Read parameterfile oon all processes for simplicity
//============================================================================
void read_parameterfile(char * fname) {

#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  char buf[500], buf1[500], buf2[500], buf3[500];
  int i,j,nt;
  int id[MAXTAGS];
  int errorFlag = 0;

  nt = 0;

  //=======================================
  // Get particles from file instead of
  // generating them in the code
  //=======================================

#ifdef READICFROMFILE
  strcpy(tag[nt], "ReadParticlesFromFile");
  addr[nt] = &ReadParticlesFromFile;
  id[nt++] = INT;

  strcpy(tag[nt], "InputParticleFileDir");
  addr[nt] = InputParticleFileDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "InputParticleFilePrefix");
  addr[nt] = InputParticleFilePrefix;
  id[nt++] = STRING;

  strcpy(tag[nt], "TypeInputParticleFiles");
  addr[nt] = &TypeInputParticleFiles;
  id[nt++] = INT;

  strcpy(tag[nt], "NumInputParticleFiles");
  addr[nt] = &NumInputParticleFiles;
  id[nt++] = INT;

  strcpy(tag[nt], "RamsesOutputNumber");
  addr[nt] = &RamsesOutputNumber;
  id[nt++] = INT;
#else
  ReadParticlesFromFile = 0;
#endif

  //===================================
  // Read modified gravity parameters
  //===================================

  read_mg_parameters(addr, tag, id, &nt);

  //===================================
  // Read the standard parameters
  //===================================

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileBase");
  addr[nt] = FileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputRedshiftFile");
  addr[nt] = &OutputRedshiftFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "NumFilesWrittenInParallel");
  addr[nt] = &NumFilesWrittenInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "UseCOLA");
  addr[nt] = &UseCOLA;
  id[nt++] = INT;

  strcpy(tag[nt], "Buffer");
  addr[nt] = &Buffer;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nmesh");
  addr[nt] = &Nmesh;
  id[nt++] = INT;

  strcpy(tag[nt], "Nsample");
  addr[nt] = &Nsample;
  id[nt++] = INT;

  strcpy(tag[nt], "Box");
  addr[nt] = &Box;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Init_Redshift");
  addr[nt] = &Init_Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Seed");
  addr[nt] = &Seed;
  id[nt++] = INT;

  strcpy(tag[nt], "SphereMode");
  addr[nt] = &SphereMode;
  id[nt++] = INT;

  strcpy(tag[nt], "WhichSpectrum");
  addr[nt] = &WhichSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "WhichTransfer");
  addr[nt] = &WhichTransfer;
  id[nt++] = INT;

  strcpy(tag[nt], "FileWithInputSpectrum");
  addr[nt] = FileWithInputSpectrum;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputTransfer");
  addr[nt] = FileWithInputTransfer;
  id[nt++] = STRING;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaBaryon");
  addr[nt] = &OmegaBaryon;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "HubbleParam");
  addr[nt] = &HubbleParam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "PrimordialIndex");
  addr[nt] = &PrimordialIndex;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
  addr[nt] = &InputSpectrum_UnitLength_in_cm;
  id[nt++] = FLOAT;

#ifndef GAUSSIAN
  strcpy(tag[nt], "Fnl_Redshift");
  addr[nt] = &Fnl_Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Fnl");
  addr[nt] = &Fnl;
  id[nt++] = FLOAT;
#endif

#ifdef GENERIC_FNL
  strcpy(tag[nt], "FileWithInputKernel");
  addr[nt] = FileWithInputKernel;
  id[nt++] = STRING;
#endif

#ifdef LIGHTCONE
  strcpy(tag[nt], "Origin_x");
  addr[nt] = &Origin_x;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Origin_y");
  addr[nt] = &Origin_y;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Origin_z");
  addr[nt] = &Origin_z;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nrep_neg_x");
  addr[nt] = &Nrep_neg_x;
  id[nt++] = INT;

  strcpy(tag[nt], "Nrep_pos_x");
  addr[nt] = &Nrep_pos_x;
  id[nt++] = INT;

  strcpy(tag[nt], "Nrep_neg_y");
  addr[nt] = &Nrep_neg_y;
  id[nt++] = INT;

  strcpy(tag[nt], "Nrep_pos_y");
  addr[nt] = &Nrep_pos_y;
  id[nt++] = INT;

  strcpy(tag[nt], "Nrep_neg_z");
  addr[nt] = &Nrep_neg_z;
  id[nt++] = INT;

  strcpy(tag[nt], "Nrep_pos_z");
  addr[nt] = &Nrep_pos_z;
  id[nt++] = INT;
#endif

  if((fd = fopen(fname, "r"))) {
    fflush(stdout);
    while(!feof(fd)) {
      buf[0] = 0;
      fgets(buf, 500, fd);

      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
      if(buf1[0] == '%') continue;

      for(i = 0, j = -1; i < nt; i++) {
        if(strcmp(buf1, tag[i]) == 0)  {
          j = i;
          tag[i][0] = 0;
          break;
        }
      }

      if(j >= 0) {
        switch (id[j]) {
          case FLOAT:
            *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
            strcpy((char *)addr[j], buf2);
            break;
          case INT:
            *((int *) addr[j]) = atoi(buf2);
            break;
        }
      } else {
        if(ThisTask == 0) fprintf(stdout,"\nERROR: In file %s:  Tag '%s' not allowed or multiple defined.\n\n",fname,buf1);
        errorFlag = 1;
      }
    }
    fclose(fd);
    for(i = 0; i < nt; i++) {
      if(*tag[i]) {
        if(ThisTask == 0) fprintf(stdout, "\nERROR: I miss a value for tag '%s' in parameter file '%s'.\n\n", tag[i], fname);
        errorFlag = 1;
      }
    }
  } else {
    if(ThisTask == 0) fprintf(stdout,"\nERROR: Parameter file '%s' not found.\n\n",fname);
    FatalError((char *)"read_param.c Parameter file not found");
  }

  if (NumFilesWrittenInParallel < 1) {
    if (ThisTask == 0) {
      printf("\nERROR: `NumFilesWrittenInParallel' is less than 1 so no processors will be writing out the data.\n");
      printf("        Please set NumFileWrittenInParallel > 0'.\n\n");
    }
    FatalError((char *)"read_param.c NumFileWrittenInParallel < 1");
  }

  if (NTask < NumFilesWrittenInParallel) {
    if (ThisTask == 0) {
      printf("\nWARNING: Number of processors smaller than `NumFilesWrittenInParallel'.\n");
      printf("         Setting NumFileWrittenInParallel = Number of processors'.\n\n");
    }
  }

  // Check the run parameters to ensure compatible gaussian/non-gaussian options
  if((WhichSpectrum != 0) && (WhichTransfer !=0)) {
    if (ThisTask == 0) {
      printf("\nERROR: You are running with both power spectrum and transfer function.\n");
      printf("       Please select the appropriate one.\n\n");
    }
    FatalError((char *)"read_param.c can't run with both power spectrum and transfer function");
  }

#ifdef GAUSSIAN
  if((WhichSpectrum == 0) && (WhichTransfer == 0)) {
    if (ThisTask == 0) {
      printf("\nERROR: You are running with neither power spectrum nor transfer function.\n");
      printf("       Please set either WhichSpectrum or WhichTransfer to a non-zero value.\n\n");
    }
    FatalError((char *)"read_param.c can't run with neither power spectrum nor transfer function");
  }
#else 
  if((WhichSpectrum != 0) || (WhichTransfer == 0)) { 
    if (ThisTask == 0) {
      printf("\nERROR: Non-Gaussian models require the transfer function as input.\n");
      printf("       Switch WhichSpectrum to zero and select a type of transfer function in the input parameter file.\n\n");
    } 
    FatalError((char *)"read_param.c Non-Gaussian models require the transfer function");
  }
#endif
#ifdef LOCAL_FNL 
  if(PrimordialIndex != 1.0) {
    if (ThisTask == 0) printf("\nERROR: Local non-gaussianity with tilted power spectrum requires the GENERIC_FNL option in the Makefile\n\n");
    FatalError((char *)"read_param.c Local non-gaussianity with tilted power spectrum requires the GENERIC_FNL define");
  }
#endif 
#ifdef ORTHO_FNL 
  if(PrimordialIndex != 1.0) {
    if (ThisTask == 0) printf("\nERROR: Orthogonal non-gaussianity with tilted power spectrum requires the GENERIC_FNL option in the Makefile\n\n"); 
    FatalError((char *)"read_param.c  Orthogonal non-gaussianity with tilted power spectrum requires the GENERIC_FNL define");
  }
#endif 
#ifdef EQUIL_FNL 
  if(PrimordialIndex != 1.0) {
    if (ThisTask == 0) printf("\nERROR: Equilateral non-gaussianity with tilted power spectrum requires the GENERIC_FNL option in the Makefile\n\n");
    FatalError((char *)"read_param.c Equilateral non-gaussianity with tilted power spectrum requires the GENERIC_FNL define");
  }
#endif 

#ifndef GAUSSIAN
  if (Fnl_Redshift < Init_Redshift) {
    if (ThisTask == 0) printf("\nERROR: Fnl_Redshift must be >= Initial Redshift for us to generate non-Gaussian initial conditions.\n\n");
    FatalError((char *)"read_param.c Fnl_Redshift must be >= Initial Redshift for non-gaussian IC");
  }
#endif

  if(errorFlag) {
    MPI_Finalize();
    exit(1);
  }

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS

  return;
}
