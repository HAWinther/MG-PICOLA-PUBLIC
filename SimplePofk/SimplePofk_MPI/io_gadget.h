#ifndef IO_GADGET
#define IO_GADGET

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void read_gadget_float_vector(FILE *fd, void *buffer, int dim, ptrdiff_t npart);
void read_gadget_int_vector(FILE *fd, void *buffer, int dim, ptrdiff_t npart);
void read_gadget_header(FILE *fd);
void process_particle_data(double *pos, double *vel, ptrdiff_t npart_now, double boxsize, struct Particles *part, int Nmesh);

//=======================================================
// Read binary method
//=======================================================

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream){
  size_t nread;
  if((nread = fread(ptr, size, nmemb, stream)) != nmemb){
    printf("I/O error (fread) nread = %li size = %li nmemb = %li\n", nread, size, nmemb);
    exit(1);
  }
  return nread;
}

//=======================================================
// Gadget header
//=======================================================

struct io_header{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char fill[60];
} header;

//=======================================================
// Print contents of gadget header
//=======================================================

void print_header(){
  printf("\n======================================\n");
  printf("Gadget Header:\n");
  printf("======================================\n");
  printf("N = %i %i %i %i %i %i\n",header.npart[0],header.npart[1],header.npart[2],header.npart[3],header.npart[4],header.npart[5]);
  printf("M = %e %e %e %e %e %e\n",header.mass[0],header.mass[1],header.mass[2],header.mass[3],header.mass[4],header.mass[5]);
  printf("a               = %f\n",header.time);
  printf("z               = %f\n",header.redshift);
  printf("BoxSize (Mpc/h) = %f\n",header.BoxSize/1000.0);
  printf("OmegaM          = %f\n",header.Omega0);
  printf("OmegaL          = %f\n",header.OmegaLambda);
  printf("HubbleParam     = %f\n",header.HubbleParam);
  printf("======================================\n\n");
}

//=======================================================
// Read gadget header
//=======================================================

void read_gadget_header(FILE *fd){
  int tmp;
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(&header, sizeof(header), 1, fd);
  my_fread(&tmp,sizeof(int),1,fd);
}

//=======================================================
// Read gadget vector of floats of size npart * dim
//=======================================================

void read_gadget_float_vector(FILE *fd, void *buffer, int dim, ptrdiff_t npart){
  int tmp;
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(buffer, dim*sizeof(float), npart, fd);
  my_fread(&tmp,sizeof(int),1,fd);
}

//=======================================================
// Read gadget vector of ints of size npart * dim
//=======================================================

void read_gadget_int_vector(FILE *fd, void *buffer, int dim, ptrdiff_t npart){
  int tmp;
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(buffer, dim*sizeof(int), npart, fd);
  my_fread(&tmp,sizeof(int),1,fd);
}

void read_and_bin_particles_gadget(char *fileprefix, int filenum, ptrdiff_t *npart_tot, float *readbuffer, ptrdiff_t *nbuffer, struct Particles *p, int Nmesh){
  char filename[100];
  sprintf(filename, "%s%i", fileprefix, filenum);
  FILE *fp;
  int npart_now;

  // Open file
  if(mpi.ThisTask == 0) printf("Opening file %s\n", filename);
  fp = fopen(filename, "r");

  // Read header and print it for the first file only
  if(filenum == 0){
    read_gadget_header(fp);
    if(mpi.ThisTask == 0) print_header();
    *npart_tot = 0;
  } else {
    read_gadget_header(fp);
  }

  npart_now = header.npart[1];

  // Check that buffer is large enough
  if( (ptrdiff_t)(3 * npart_now) > *nbuffer){
    printf("Error: increase size of buffer. Trying to read %i > nbuffer = %td\n", 3*npart_now, *nbuffer);
    exit(1); 
  }
   
  // Verbose
  if(mpi.ThisTask == 0){
    printf("=> Reading file: %s\n", filename);
    printf("=> Filenum = %i  Npart_loc = %i  Npart_tot = %lli \n", filenum, npart_now, (long long)(*npart_tot+npart_now));
  }

  // Read positions
  read_gadget_float_vector(fp, readbuffer, 3, npart_now);

#ifdef _READVELOCITY
  // Read velocities
  float *readvelbuffer = &readbuffer[3*npart_now];
  read_gadget_float_vector(fp, readvelbuffer, 3, npart_now);
  process_particle_data((double *)readbuffer, (double *)readvelbuffer, (long long)(npart_now), (double)(header.BoxSize), p, Nmesh);
#else
  process_particle_data((double *)readbuffer, NULL, (long long)(npart_now), (double)(header.BoxSize), p, Nmesh);
#endif
  

  // Update how many particles in total we have read
  *npart_tot += (long long)(npart_now);

  fclose(fp);
}

#endif
