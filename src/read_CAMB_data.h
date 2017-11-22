
//==================================
// Format of transfer file
// Standard values below is for CAMB
//==================================
const int n_transfer_header_lines = 1;    // Number of header lines
const int ncol_transfer_file      = 13;   // Columns in file
const int transfer_col_k          = 0;    // col number of k
const int transfer_col_cdm        = 1;    // col number of T_CDM
const int transfer_col_baryon     = 2;    // col number of T_b
const int transfer_col_nu         = 5;    // col number of T_massive_nu
const int transfer_col_total      = 6;    // col number of T_total
const int NMAX_ROWS_TRANSFER      = 5000; // Maximum number of rows
//==================================

//==================================
// Format of P(k) file
// Standard values below is for CAMB
//==================================
const int n_pofk_header_lines     = 1;    // Number of header lines
const int ncol_pofk_file          = 2;    // Columns in file
const int pofk_col_k              = 0;    // col number of k
const int pofk_col_pofk           = 1;    // col number of P(k)
const int NMAX_ROWS_POFK          = 5000; // Maximum number of rows
//==================================

GSL_2D_Spline nu_transfer_function_spline;
GSL_2D_Spline cdm_transfer_function_spline;
GSL_2D_Spline baryon_transfer_function_spline;
GSL_2D_Spline total_transfer_function_spline;
GSL_Spline    total_pofk_spline;

void free_CAMB_splines(){
  Free_GSL_2D_Spline(&nu_transfer_function_spline);
  Free_GSL_2D_Spline(&cdm_transfer_function_spline);
  Free_GSL_2D_Spline(&baryon_transfer_function_spline);
  Free_GSL_2D_Spline(&total_transfer_function_spline);
  Free_GSL_Spline(&total_pofk_spline);
}

void read_and_spline_transfer_functions(char *transferfileinfofilename){

  // Temporary memory to contain the data
  double *transfer_function_nu_tmp     = my_malloc(sizeof(double)*NMAX_ROWS_TRANSFER);
  double *transfer_function_baryon_tmp = my_malloc(sizeof(double)*NMAX_ROWS_TRANSFER);
  double *transfer_function_cdm_tmp    = my_malloc(sizeof(double)*NMAX_ROWS_TRANSFER);
  double *transfer_function_total_tmp  = my_malloc(sizeof(double)*NMAX_ROWS_TRANSFER);
  double *logk_tmp                     = my_malloc(sizeof(double)*NMAX_ROWS_TRANSFER);

  // Open fileinfo file
  char filepath[400];
  int nredshift;
  FILE *fp = fopen(transferfileinfofilename,"r");
  if(fp == NULL){
    printf("Error read_and_spline_transfer_functions: cannot read [%s]\n", transferfileinfofilename);
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }

  // Read fileprefix and nfiles
  fscanf(fp, "%s %i\n", filepath, &nredshift);
  if(ThisTask == 0){
    printf("Reading transfer functions | Filedir [%s] | Reading [%i] redshift files\n", filepath, nredshift);
  }

  double *transfer_function_nu = NULL, *transfer_function_cdm = NULL, *logk = NULL,
         *transfer_function_baryon = NULL, *transfer_function_total = NULL;
  double *redshifts  = my_malloc(sizeof(double) * nredshift);
  int nk = 0;

  // Read all files
  for(int i = 0; i < nredshift; i++){
    double znow;
    int npts;

    // Read filename
    char filename[400];
    fscanf(fp,"%s %lf\n", filename, &znow);
    redshifts[i] = znow;

    // Make filename and open file. Assumes all files have the same length
    char fullfilename[400];
    sprintf(fullfilename,"%s/%s", filepath, filename);
    read_single_transfer_function_file(fullfilename, logk_tmp, transfer_function_nu_tmp,
        transfer_function_cdm_tmp, transfer_function_baryon_tmp, transfer_function_total_tmp, &npts);

    // Allocate memory
    if(i == 0){
      transfer_function_nu     = my_malloc(sizeof(double) * nredshift * npts);
      transfer_function_cdm    = my_malloc(sizeof(double) * nredshift * npts);
      transfer_function_baryon = my_malloc(sizeof(double) * nredshift * npts);
      transfer_function_total  = my_malloc(sizeof(double) * nredshift * npts);
      logk                     = my_malloc(sizeof(double) * nredshift * npts);
      nk = npts;
    }

    // Check that k-array is the same in all files as this is assumed when splining below
    if(i > 0){
      double maxerr = 0.0;
      for(int j = 0; j < nk; j++){
        double err = fabs(logk_tmp[j]/logk[j]-1.0);
        if(err > maxerr) err = maxerr;
      }
      if(maxerr > 1e-3){
        printf("Error in read_and_spline_transfer_function the k-array differs in the different files. Not built-in support for this\n");
        MPI_Abort(MPI_COMM_WORLD,1);
        exit(1);
      }
    }

    // Copy over data
    for(int j = 0; j < nk; j++){
      logk[j] = logk_tmp[j];
      transfer_function_nu[i * nk + j]     = transfer_function_nu_tmp[j];
      transfer_function_cdm[i * nk + j]    = transfer_function_cdm_tmp[j];
      transfer_function_baryon[i * nk + j] = transfer_function_baryon_tmp[j];
      transfer_function_total[i * nk + j]  = transfer_function_total_tmp[j];
    }

    if(ThisTask == 0){
      printf("z = [%8.4f] | We have [%5i] k-points | Filename: [%s]\n", znow, npts, fullfilename);
    }
  }
  fclose(fp);

  // Create spline
  Create_GSL_2D_Spline(&nu_transfer_function_spline,     logk, redshifts, transfer_function_nu,     nk, nredshift);
  Create_GSL_2D_Spline(&cdm_transfer_function_spline,    logk, redshifts, transfer_function_cdm,    nk, nredshift);
  Create_GSL_2D_Spline(&baryon_transfer_function_spline, logk, redshifts, transfer_function_baryon, nk, nredshift);
  Create_GSL_2D_Spline(&total_transfer_function_spline,  logk, redshifts, transfer_function_total,  nk, nredshift);

  // Test spline
  if(ThisTask == 0){
    printf("\nSample values neutrino transfer function:\n");
    for(int i = 0; i < nredshift; i++){
      for(int j = 0; j < nk; j++){
        if(rand() % 200 == 0){
          printf("z: %5.2f k: %8.4f   Tnu/Tnu0   / D: %10.5f \n", redshifts[i], exp(logk[j]),
              transfer_function_nu[i * nk + j]  / (transfer_function_nu[j] + 1e-20)  / growth_D_LCDMFit(1.0/(1.0+redshifts[i])));
        }
      }
    }
    printf("\nSample values CDM transfer function:\n");
    for(int i = 0; i < nredshift; i++){
      for(int j = 0; j < nk; j++){
        if(rand() % 200 == 0){
          printf("z: %5.2f k: %8.4f   Tcdm/Tcdm0 / D: %10.5f \n", redshifts[i], exp(logk[j]),
              transfer_function_cdm[i * nk + j] / transfer_function_cdm[j] / growth_D_LCDMFit(1.0/(1.0+redshifts[i])));
        }
      }
    }
    printf("==============================================================\n\n");
  }

  // Free memory
  my_free(redshifts);
  my_free(logk);
  my_free(transfer_function_nu);
  my_free(transfer_function_cdm);
  my_free(transfer_function_baryon);
  my_free(transfer_function_total);
  my_free(logk_tmp);
  my_free(transfer_function_nu_tmp);
  my_free(transfer_function_cdm_tmp);
  my_free(transfer_function_baryon_tmp);
  my_free(transfer_function_total_tmp);
}

//=============================================
// Massive nu transfer function
//=============================================
double get_nu_transfer_function(double k, double a){
  double z = 1.0/a - 1.0;
  return Lookup_GSL_2D_Spline(&nu_transfer_function_spline, log(k), z);
}

//=============================================
// CDM transfer function
//=============================================
double get_cdm_transfer_function(double k, double a){
  double z = 1.0/a - 1.0;
  return Lookup_GSL_2D_Spline(&cdm_transfer_function_spline, log(k), z);
}

//=============================================
// Baryon transfer function
//=============================================
double get_baryon_transfer_function(double k, double a){
  double z = 1.0/a - 1.0;
  return Lookup_GSL_2D_Spline(&baryon_transfer_function_spline, log(k), z);
}

//=============================================
// Weighted CDM+Baryon transfer function
// Since baryons are CDM in the simulations
// this is the one to use
//=============================================
double get_cdm_baryon_transfer_function(double k, double a){
  return (OmegaBaryon * get_baryon_transfer_function(k,a) + OmegaCDM * get_cdm_transfer_function(k,a))/(OmegaBaryon + OmegaCDM);
}

//=============================================
// Total (CDM+b+nu) transfer function
//=============================================
double get_total_transfer_function(double k, double a){
  double z = 1.0/a - 1.0;
  return Lookup_GSL_2D_Spline(&total_transfer_function_spline, log(k), z);
}

//=============================================
// Total (CDM+b+nu) power-spectrum 
// As written this is P(k) / A_s
// where A_s is the primordial amplitude
//=============================================
double get_total_power_spectrum(double k){
  return 2.0 * M_PI * M_PI * pow(get_total_transfer_function(k, 1.0), 2) * (k * HubbleParam) 
   * pow(k * HubbleParam / 0.05, PrimordialIndex - 1.0) * pow(HubbleParam, 3);
}

//=============================================
// (CDM+b) power-spectrum 
//=============================================
double get_cdm_baryon_power_spectrum(double k){
  return get_total_power_spectrum(k) * pow( get_cdm_baryon_transfer_function(k, 1.0) / get_total_transfer_function(k, 1.0), 2 );
}

//=============================================
// Neutrino power-spectrum 
//=============================================
double get_neutrino_power_spectrum(double k){
  return get_total_power_spectrum(k) * pow( get_nu_transfer_function(k, 1.0) / get_total_transfer_function(k, 1.0), 2 );
}

//=============================================
// Read a single (CAMB) transfer function file
// and return CDM and massive-nu transfer functions
//=============================================
void read_single_transfer_function_file(char *filename, double *logk,
    double *transfer_function_nu, double *transfer_function_cdm, double *transfer_function_baryon, double *transfer_function_total, int *npts){
  double tmp[ncol_transfer_file];

  // Open CAMB transfer function file for reading
  FILE *fp = fopen(filename, "r");
  if(fp == NULL){
    printf("Error ed_camb_transfer_function_file: cannot open [%s]\n", filename);
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  // Read first line
  for(int i = 0; i < n_transfer_header_lines; i++){
    size_t linelen = 0;
    char *firstline = NULL;
    getline(&firstline, &linelen, fp);
    free(firstline);
  }

  // Read data row by row
  int nrows = 0;
  for(;;){
    if( fscanf(fp,"%lf", &tmp[0]) == EOF ) break;
    for(int i = 1; i < ncol_transfer_file; i++)
      fscanf(fp, "%lf", &tmp[i]);

    // 0: k/h   1: CDM      2: baryon   3: photon  4: nu     5: mass_nu  6: total 
    // 7: no_nu 8: total_de 9: Weyl    10:        11: v_CDM 12: v_b     13: v_b-v_c 
    logk[nrows]                     = log(tmp[transfer_col_k]);
    transfer_function_cdm[nrows]    = tmp[transfer_col_cdm];
    transfer_function_baryon[nrows] = tmp[transfer_col_baryon];
    transfer_function_nu[nrows]     = tmp[transfer_col_nu];
    transfer_function_total[nrows]  = tmp[transfer_col_total];
    nrows++;
  }
  fclose(fp);

  *npts = nrows;
}

