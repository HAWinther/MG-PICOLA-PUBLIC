#ifndef IO_RAMSES
#define IO_RAMSES

void process_particle_data(double *pos, int npart, double boxsize);

//======================================
// Read method for RAMSES
//======================================

int read_int(FILE* fp){
  int tmp, skip;
  fread(&skip, sizeof(int), 1, fp);
  fread(&tmp, sizeof(int), 1, fp);
  fread(&skip, sizeof(int), 1, fp);
  return tmp;
}

void read_int_vec(FILE* fp, int *buffer, int n){
  int skip;
  fread(&skip, sizeof(int), 1, fp);
  fread(buffer, sizeof(int), n, fp);
  fread(&skip, sizeof(int), 1, fp);
}

void read_double_vec(FILE* fp, double *buffer, int n){
  int skip;
  fread(&skip, sizeof(int), 1, fp);
  fread(buffer, sizeof(double), n, fp);
  fread(&skip, sizeof(int), 1, fp);
}

std::string int_to_ramses_string(int i){
  char *filenameout = new char[200];
  std::string filename;
  if(i<10)
    sprintf(filenameout,"0000%i",i);
  else if(i<100)
    sprintf(filenameout,"000%i",i);
  else if(i<1000)
    sprintf(filenameout,"00%i",i);
  else if(i<10000)
    sprintf(filenameout,"0%i",i);
  else{
    std::cout << "add_ramses_ending not implemented for 10000 <= " << i << std::endl;
    exit(1);
  }
  filename = std::string(filenameout);
  delete[] filenameout;
  return filename;
}

//======================================
// Read method for RAMSES
//======================================

void read_and_bin_particles_ramses(std::string filebase, int filenum, int *npart_tot, double *readbuffer, int *nbuffer){
  FILE *fp;
  std::string filename;

  //======================================
  // Open file
  //======================================
  filename = filebase + int_to_ramses_string(filenum); 
  std::cout << "Opening file: " << filename << std::endl;
  fp = fopen(filename.c_str(), "r");

  //======================================
  // Read header
  //======================================
  int ncpu_loc, ndim_loc, npart_loc, localseed_loc[4], nstar_tot_loc, mstar_tot_loc[2], mstar_lost_loc[2], nsink_loc;
  ncpu_loc      = read_int(fp);
  ndim_loc      = read_int(fp);
  npart_loc     = read_int(fp);
  read_int_vec(fp, localseed_loc , 4);
  nstar_tot_loc = read_int(fp);
  read_int_vec(fp, mstar_tot_loc , 2);
  read_int_vec(fp, mstar_lost_loc, 2);
  nsink_loc     = read_int(fp);
  (void) ncpu_loc;
  (void) nsink_loc;
  (void) nstar_tot_loc;

  // Sanity check
  if(ndim_loc*npart_loc > *nbuffer){
    std::cout << "Error: increase size of buffer. Trying to read " << 3*npart_loc << " > nbuffer = " << *nbuffer << std::endl;
    exit(1); 
  }

  // Verbose
  std::cout << "=> Reading file: " << filename << std::endl;
  std::cout << "=> Filenum = " << filenum << " Npart_loc = " << npart_loc << " Npart_tot = " << *npart_tot+npart_loc << std::endl;
  *npart_tot += npart_loc;

  //======================================
  // Read particles
  //======================================
  read_double_vec(fp, &readbuffer[0], npart_loc);
  read_double_vec(fp, &readbuffer[npart_loc], npart_loc);
  read_double_vec(fp, &readbuffer[2*npart_loc], npart_loc);

  process_particle_data(readbuffer, npart_loc, 0.0);

  fclose(fp);
}

#endif
