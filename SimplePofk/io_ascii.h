#ifndef IO_ASCII
#define IO_ASCII

void process_particle_data(double *pos, int npart, double boxsize);

//======================================
// Read method for ascii
//======================================

void read_and_bin_particles_ascii(std::string filebase, int filenum, int *npart_tot, double *readbuffer, int *nbuffer){
  FILE *fp;
  std::string filename;

  //======================================
  // Open file 
  //======================================
  filename = filebase; 
  std::cout << "Opening file: " << filename << std::endl;
  fp = fopen(filename.c_str(),"r");

  // Read npart
  int npart_loc;
  fscanf(fp,"%i\n",&npart_loc);
  *npart_tot = npart_loc;

  // Check if buffer is large enough
  if(3*npart_loc > *nbuffer){
    std::cout << "Error: increase size of buffer. Trying to read " << 3*npart_loc << " > nbuffer = " << *nbuffer << std::endl;
    exit(1);
  }

  // Verbose
  std::cout << "=> Reading file: " << filename << " npart: " << npart_loc << std::endl;

  //======================================
  // Read particles
  //======================================

  for(int i = 0;i < npart_loc; i++){
    double tmp;
    //fscanf(fp,"%lf %lf %lf %lf %lf %lf\n", &readbuffer[3*i], &readbuffer[3*i+1], &readbuffer[3*i+2], &tmp, &tmp, &tmp);
    fscanf(fp,"%lf %lf %lf %lf\n", &readbuffer[3*i], &readbuffer[3*i+1], &readbuffer[3*i+2], &tmp);
  }

  // Calculate boxsize as max x,y,z
  double boxsize = 0.0;
  for(int i=0;i<3*npart_loc;i++)
    if(readbuffer[i] > boxsize) boxsize = readbuffer[i];
  std::cout << "=> Boxsize: " << boxsize << std::endl;

  process_particle_data(readbuffer, npart_loc, boxsize);

  fclose(fp);
}

#endif
