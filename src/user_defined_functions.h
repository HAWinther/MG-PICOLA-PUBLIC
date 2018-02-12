#include "msg.h"
#include "timer.h"

//=========================================================//
//                                                         //
// MG-PICOLA written by Hans Winther (ICG Portsmouth) 2017 //
//                                                         //
// If implementing a new model the functions one need to   //
// change in this file is:                                 //
//                                                         //
// - hubble(a) dhubbleda(a) (if modified cosmology)        //
//                                                         //
// - GeffofG(a,k) aka mu(a,k) the mod to the Poisson eq.   //
//                                                         //
// - The second order equivalents: second_order_kernel     //
//   or Factor_2LPT (if only time-dependent 2LPT)          //
//                                                         //
// - Add new parameters to vars.h / vars.c and read them   //
//   from file in the read_param routine below             //
//                                                         //
// - Screening method is choosen in ComputeFifthForce      //
//   and the screening function is set in the function     //
//   screening_factor_X where X is density, potential      //
//   or gradient.                                          //
//                                                         //
// - The rest of the code should not need any changes      //
//                                                         //
// To use f(R) gravity with scale-dependent growth use     //
// -DFOFRGRAVITY and -DSCALEDEPENDENT                      //
//                                                         //
// To use DGP gravity use -DDGPGRAVITY                     //
//                                                         //
//=========================================================//

//=========================================================//
// This functions is called right after parameters have    //
// been read                                               //
// This is run before starting to compute growth-factors   //                           
//=========================================================//
void init_modified_version(){

  //==========================================================
  // This flag must be set if any of the screening methods
  // are to be used and/or if the MG potential must be computed
  // For models where phi = C(a) Phi_N, i.e. Geff(a) = 1 + C(a)
  // then we don't need this
  //==========================================================
  allocate_mg_arrays = 1;

#ifdef MASSIVE_NEUTRINOS

#ifndef SCALEDEPENDENT
  printf("Massive neutrinos must be compiled with SCALEDEPENDENT\n");
  MPI_Abort(MPI_COMM_WORLD,1);
  exit(1);
#endif

  // Set neutrino density parameter
  OmegaNu  = nu_SumMassNuEV / (93.14 * HubbleParam * HubbleParam);
  OmegaCDM = Omega - OmegaBaryon - OmegaNu;

  // For very small values it does not contribute so turn it off
  if(OmegaNu < 1e-6) {
    OmegaNu  = 0.0;
    OmegaCDM = Omega - OmegaBaryon;
    nu_include_massive_neutrinos = 0;
  }

  if(nu_include_massive_neutrinos){
    if(ThisTask == 0){
      printf("==============================================================\n");
      printf("Running with massive neutrinos SumMass = %f eV  OmegaNu = %f\n", nu_SumMassNuEV, OmegaNu);
      printf("==============================================================\n");
    }
  }

#endif

#if defined(FOFRGRAVITY)

  if(ThisTask == 0 && modified_gravity_active){
    printf("============================================\n");
    printf("Running with Modified Gravity, f(R) gravity \n");
    printf("============================================\n");
    printf("  f(R0)                   = %f\n", fofr0);
    printf("  n                       = %f\n", nfofr);
    printf("  Modified gravity active = %i\n", modified_gravity_active);
    printf("  Include screening       = %i\n", include_screening);
    printf("  Use LCDM growthfactor   = %i\n", use_lcdm_growth_factors);
    printf("  Input P(k) is for LCDM  = %i\n", input_pofk_is_for_lcdm);
    printf("  Sigma8 is for LCDM      = %i\n", input_sigma8_is_for_lcdm);
    printf("\n");
    fflush(stdout);
  }

#elif defined(DGPGRAVITY)

  if(ThisTask == 0 && modified_gravity_active){
    printf("============================================\n");
    printf("Running with Modified Gravity, nDGP gravity \n");
    printf("============================================\n");
    printf("  rcH0                    = %f\n", rcH0_DGP);
    printf("  Rsmooth (Mpc/h)         = %f\n", Rsmooth_global);
    printf("  Modified gravity active = %i\n", modified_gravity_active);
    printf("  Include screening       = %i\n", include_screening);
    printf("  Use LCDM growthfactor   = %i\n", use_lcdm_growth_factors);
    printf("  Input P(k) is for LCDM  = %i\n", input_pofk_is_for_lcdm);
    printf("  Sigma8 is for LCDM      = %i\n", input_sigma8_is_for_lcdm);
    printf("\n");
    fflush(stdout);
  }

#elif defined(MBETAMODEL)

  if(ThisTask == 0 && modified_gravity_active){
    printf("============================================\n");
    printf("Running with Modified Gravity, (m(a),b(a))  \n");
    printf("The example model used is the symmetron     \n");
    printf("============================================\n");
    printf("  Modified gravity active = %i\n", modified_gravity_active);
    printf("  Include screening       = %i\n", include_screening);
    printf("  Use LCDM growthfactor   = %i\n", use_lcdm_growth_factors);
    printf("  Input P(k) is for LCDM  = %i\n", input_pofk_is_for_lcdm);
    printf("  Sigma8 is for LCDM      = %i\n", input_sigma8_is_for_lcdm);
    printf("  assb_symm               = %f\n", assb_symm);
    printf("  beta_symm               = %f\n", beta_symm);
    printf("  range_symm              = %f Mpc/h\n", range_symm);
    printf("\n");
    fflush(stdout);
  }

  // For (m,beta) models we need to compute Phi_crit(a) and store in spline
  if(modified_gravity_active) compute_phi_of_a();

#elif defined(KMOFLAGE)

  // Model not implemented...

#elif defined(BRANSDICKE)

  // We don't need any MG arrays here
  allocate_mg_arrays = 0;

  if(ThisTask == 0){
    printf("============================================\n");
    printf("Running with the Jordan-Brans-Dicke model   \n");
    printf("============================================\n");
    printf("  Modified gravity active = %i\n", modified_gravity_active);
    printf("  Input P(k) is for LCDM  = %i\n", input_pofk_is_for_lcdm);
    printf("  Sigma8 is for LCDM      = %i\n", input_sigma8_is_for_lcdm);
    printf("\n");
    fflush(stdout);
  }

  // Solve background and get the Hubble parameter today
  double h;
  JBD_Solve_Background(wBD, Omegah2, Omegavh2, Omegarh2, &h);
  
  // Set the values we find
  HubbleParam = h;
  Omega       = Omegah2 / HubbleParam / HubbleParam;

  if(ThisTask == 0){
    printf("Calling JBD Background Solver: \n");
    printf(" OmegaMh2 = %f\n", Omegah2);
    printf(" OmegaRh2 = %f\n", Omegarh2);
    printf(" OmegaVh2 = %f\n", Omegavh2);
    printf(" wBD      = %f\n", wBD);
    printf(" We find OmegaM     : %f\n", Omega);
    printf(" We find HubbleParam: %f\n", HubbleParam);
    printf("\nTest of JBD splines: \n");
    printf("GeffG(a=1.0) = %7.3f\n", GeffoverG(1.0, 0.0));
    printf("GeffG(a=0.5) = %7.3f\n", GeffoverG(0.5, 0.0));
    printf("GeffG(a=0.1) = %7.3f\n", GeffoverG(0.1, 0.0));
    printf("H(a=1.0) = %7.3f  HLCDM(a=1.0) = %7.3f\n", hubble(1.0), 1.0);
    printf("H(a=0.5) = %7.3f  HLCDM(a=0.5) = %7.3f\n", hubble(0.5), sqrt(Omega/(0.5*0.5*0.5) + 1.0 - Omega));
    printf("H(a=0.1) = %7.3f  HLCDM(a=0.1) = %7.3f\n", hubble(0.1), sqrt(Omega/(0.1*0.1*0.1) + 1.0 - Omega));
    fflush(stdout);
    fflush(stdout);
  }

#endif

#if defined(EQUATIONOFSTATE_PARAMETRIZATION)
  
  //=================================================================
  // If this is active then we use the (w0,wa) background
  // Not compatible with JBD and other models where we don't have a
  // LCDM background
  //=================================================================
  if(ThisTask == 0) {
    printf("\n");
    printf("===========================\n");
    printf("Using (w0,wa) for background  \n");
    printf("===========================\n");
    printf("w_0   = %f\n", w_0);
    printf("w_a   = %f\n", w_a);
    printf("\n");
    printf("H(a=0.5) = %7.3f  HLCDM(a=0.5) = %7.3f\n", hubble(0.5), sqrt(Omega/pow(0.5,3) + 1 - Omega));
    printf("H(a=0.1) = %7.3f  HLCDM(a=0.1) = %7.3f\n", hubble(0.1), sqrt(Omega/pow(0.1,3) + 1 - Omega));
    printf("dHda(a=0.5) = %7.3f  dHLCDMda(a=0.5) = %7.3f\n", 
        dhubbleda(0.5), 1.0/(2.0 * sqrt(Omega/pow(0.5,3) + 1 - Omega)) * ( -3.0 * Omega / pow(0.5,4) ));
    printf("dHda(a=0.1) = %7.3f  dHLCDMda(a=0.1) = %7.3f\n", 
        dhubbleda(0.1), 1.0/(2.0 * sqrt(Omega/pow(0.1,3) + 1 - Omega)) * ( -3.0 * Omega / pow(0.1,4) ));
    fflush(stdout);
  }
#endif

  //=================================================================
  // Compute first and second order growth-factors and spline them up
  //=================================================================
  if(ThisTask == 0) {
    printf("\n");
    printf("===========================\n");
    printf("Compute growth-factors...  \n");
    printf("===========================\n");
    fflush(stdout);
  }
  solve_for_growth_factors();

}

//=========================================================//
// Read MG parameters from file                            //
// The parameters set here needs to be defined in vars.*   //
//=========================================================//
void read_mg_parameters(void **addr, char (*tag)[50], int *id, int (*nt)){

#define FLOAT  1
#define STRING 2
#define INT    3

  strcpy(tag[(*nt)], "modified_gravity_active");
  addr[(*nt)] = &modified_gravity_active;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "use_lcdm_growth_factors");
  addr[(*nt)] = &use_lcdm_growth_factors;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "input_pofk_is_for_lcdm");
  addr[(*nt)] = &input_pofk_is_for_lcdm;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "input_sigma8_is_for_lcdm");
  addr[(*nt)] = &input_sigma8_is_for_lcdm;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "include_screening");
  addr[(*nt)] = &include_screening;
  id[(*nt)++] = INT;

#ifdef MASSIVE_NEUTRINOS
  strcpy(tag[(*nt)], "nu_include_massive_neutrinos");
  addr[(*nt)] = &nu_include_massive_neutrinos;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "nu_SumMassNuEV");
  addr[(*nt)] = &nu_SumMassNuEV;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "nu_FilenameTransferInfofile");
  addr[(*nt)] = &nu_FilenameTransferInfofile;
  id[(*nt)++] = STRING;
#endif

#ifdef MATCHMAKER_HALOFINDER

  strcpy(tag[(*nt)], "mm_run_matchmaker");
  addr[(*nt)] = &mm_run_matchmaker;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "mm_output_pernode");
  addr[(*nt)] = &mm_output_pernode;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "mm_output_format");
  addr[(*nt)] = &mm_output_format;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "mm_min_npart_halo");
  addr[(*nt)] = &mm_min_npart_halo;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "mm_linking_length");
  addr[(*nt)] = &mm_linking_length;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "mm_dx_extra_mpc");
  addr[(*nt)] = &mm_dx_extra_mpc;
  id[(*nt)++] = FLOAT;

#endif

#ifdef COMPUTE_POFK

  //=======================================
  // Compute power-spectrum within code
  //=======================================

  strcpy(tag[(*nt)], "pofk_compute_every_step");
  addr[(*nt)] = &pofk_compute_every_step;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "pofk_compute_rsd_pofk");
  addr[(*nt)] = &pofk_compute_rsd_pofk;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "pofk_nbins");
  addr[(*nt)] = &pofk_nbins;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "pofk_bintype");
  addr[(*nt)] = &pofk_bintype;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "pofk_subtract_shotnoise");
  addr[(*nt)] = &pofk_subtract_shotnoise;
  id[(*nt)++] = INT;

  strcpy(tag[(*nt)], "pofk_kmin");
  addr[(*nt)] = &pofk_kmin;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "pofk_kmax");
  addr[(*nt)] = &pofk_kmax;
  id[(*nt)++] = FLOAT;

#endif

#if defined(EQUATIONOFSTATE_PARAMETRIZATION)
  strcpy(tag[(*nt)], "w_0");
  addr[(*nt)] = &w_0;
  id[(*nt)++] = FLOAT;
  
  strcpy(tag[(*nt)], "w_a");
  addr[(*nt)] = &w_a;
  id[(*nt)++] = FLOAT;
#endif

#if defined(FOFRGRAVITY) 

  strcpy(tag[(*nt)], "fofr0");
  addr[(*nt)] = &fofr0;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "nfofr");
  addr[(*nt)] = &nfofr;
  id[(*nt)++] = FLOAT;

#elif defined(MBETAMODEL)

  strcpy(tag[(*nt)], "range_symm");
  addr[(*nt)] = &range_symm;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "assb_symm");
  addr[(*nt)] = &assb_symm;
  id[(*nt)++] = FLOAT;
  
  strcpy(tag[(*nt)], "beta_symm");
  addr[(*nt)] = &beta_symm;
  id[(*nt)++] = FLOAT;

#elif defined(DGPGRAVITY)

  strcpy(tag[(*nt)], "Rsmooth");
  addr[(*nt)] = &Rsmooth_global;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "rcH0_DGP");
  addr[(*nt)] = &rcH0_DGP;
  id[(*nt)++] = FLOAT;

#elif defined(BRANSDICKE)
  
  strcpy(tag[(*nt)], "wBD");
  addr[(*nt)] = &wBD;
  id[(*nt)++] = FLOAT;
  
  strcpy(tag[(*nt)], "Omegah2");
  addr[(*nt)] = &Omegah2;
  id[(*nt)++] = FLOAT;
  
  strcpy(tag[(*nt)], "Omegar2");
  addr[(*nt)] = &Omegarh2;
  id[(*nt)++] = FLOAT;
  
  strcpy(tag[(*nt)], "Omegav2");
  addr[(*nt)] = &Omegavh2;
  id[(*nt)++] = FLOAT;

#else

  // Define your parameters here...

#endif
}

//=========================================================//
// Hubble function and derivative                          //
// NB: if these are modified then what we call LCDM growth //
// factors will also be modified accordingly.              //  
//=========================================================//
double hubble(double a){

#if defined(BRANSDICKE)

  return JBD_Hubble_of_a(a);

#elif defined(EQUATIONOFSTATE_PARAMETRIZATION)  

  // w(a) = w0 + wa(1-a) parametrization
  return sqrt( Omega/(a*a*a) + (1.0 - Omega) * exp( 3.0 * w_a * (a-1) - 3*(1 + w_0 + w_a) * log(a))  );

#else

  return sqrt(Omega / (a*a*a) + 1.0 - Omega);

#endif
}

double dhubbleda(double a){

#if defined(BRANSDICKE)

  return JBD_dHubbleda_of_a(a);

#elif defined(EQUATIONOFSTATE_PARAMETRIZATION)  

  // w(a) = w0 + wa(1-a) parametrization
  return 1.0/(2.0 * hubble(a)) * ( -3.0 * Omega / (a*a*a*a) - 3.0 * (1.0 - Omega) * (1.0 + w_0 + w_a*(1-a) ) / a );

#else

  return 1.0/(2.0 * hubble(a)) * ( -3.0 * Omega / (a*a*a*a) );

#endif
}

#ifdef DGPGRAVITY
//=========================================================//
// The DGP beta function (Geff/G(a) = 1  + 1/3beta(a))     //
//=========================================================//
double beta_DGP(double a){
  return 1.0 + 2.0 * rcH0_DGP * ( hubble(a) +  a * dhubbleda(a) / 3.0 );
}
#endif

//=========================================================//
// The beta(a) function. When inside range and unscreened  //
// the fifth force is 2beta^2 times gravity                //
//=========================================================//
double beta_of_a(double a){
#if defined(FOFRGRAVITY)

  return 1.0/sqrt(6.0);

#elif defined(MBETAMODEL)
  
  // For the symmetron model we have
  if(a > assb_symm){
    return beta_symm * sqrt(1.0 - pow3(assb_symm/a));
  } else {
    return 0.0;
  }

  // For f(R) we have
  // return 1.0/sqrt(6.0);

#else

  return 0.0;

#endif
}

//=========================================================//
// The m(a) function (acctually m^2(a)/H0^2)               //
// The range of the force is 1/m(a). Not used for DGP      //
//=========================================================//
double mass2_of_a(double a){
#if defined(FOFRGRAVITY)

  // This is m^2/H0^2
  double a3    = a * a * a;
  double fac   = Omega/a3 + 4.0 * (1.0-Omega);
  double fac0  = Omega    + 4.0 * (1.0-Omega);
  double mass2 = fac0 * pow( fac / fac0, nfofr + 2.0) / ((1.0 + nfofr) * fofr0);
  return mass2; 

#elif defined(MBETAMODEL)

  // For the symmetron we have
  if(a > assb_symm){
    return 0.5 * pow2(2998.0 / range_symm ) * (1.0 - pow3(assb_symm/a));
  } else {
    return 1.0; // Unimporant as beta == 0, but non-zero value 
                // useful to avoid errors in computing phi(a)
  }

  // For Hu-Sawicky f(R) we have
  // double a3    = a * a * a;
  // double fac   = Omega/a3 + 4.0 * (1.0-Omega);
  // double fac0  = Omega    + 4.0 * (1.0-Omega);
  // double mass2 = fac0 * pow( fac / fac0, nfofr + 2.0) / ((1.0 + nfofr) * fofr0);
  // return mass2; 

#else
  
  return 0.0;

#endif
}

//=========================================================//
// The derivative dm^2(a)/da needed for the gamma2 factor
//=========================================================//
double dmass2_of_ada(double a){
#if defined(FOFRGRAVITY)

  // This is dm^2(a)/da / H0^2
  double a3    = a * a * a;
  double fac   = Omega/a3 + 4.0 * (1.0-Omega);
  double fac0  = Omega    + 4.0 * (1.0-Omega);
  double mass2 = fac0 * (nfofr + 2.0) * pow( fac / fac0, nfofr + 1.0) / ((1.0 + nfofr) * fofr0 * fac0);
  mass2 *= -3.0*Omega/a3/a;
  return mass2; 

#elif defined(MBETAMODEL)
  
  // For the symmetron we have
  if(a > assb_symm){
    return 0.5 * pow2(2998.0 / range_symm ) * ( + 3.0 * pow3(assb_symm/a) / a);
  } else {
    return 0.0;
  }

  // For Hu-Sawicky f(R) we have
  // double a3    = a * a * a;
  // double fac   = Omega/a3 + 4.0 * (1.0-Omega);
  // double fac0  = Omega    + 4.0 * (1.0-Omega);
  // double mass2 = fac0 * (nfofr + 2.0) * pow( fac / fac0, nfofr + 1.0) / ((1.0 + nfofr) * fofr0 * fac0);
  // mass2 *= -3.0*Omega/a3/a;
  // return mass2; 

#else

  return 0.0;

#endif
}

//=========================================================//
// This is phi_crit(a) = [9 Omegam] / [2beta(a)]*Int_ai^a  // 
//  [ beta_of_a(A) / (mass2_of_a(A) * A * A * A * A) dA ]  //
// Needed for models defined by beta(a) and m(a). This is  //
// phi(a)/2beta(a) in the notation of Brax et al. This is  //
// computed by the code, but if the analytical expression  //
// is known it's better to just define it here             //
//=========================================================//
#if defined(MBETAMODEL)
double phi_of_a(double a){

  // E.g. for the symmetron we have
  // return 3.0 * Omega / pow3(assb_symm) * pow2(range_symm / 2998.0);

  // Use what we computed in src/cosmo.c::compute_phi_of_a()
  return phi_of_a_from_spline(a);
}
#endif

//=========================================================//
// The general linear MG equation is                       //
// 1/a^2 D^2 phi = m^2 phi + coupling(a) * 4piG a^2 * drho //
// where the total force is F_tot = DPhi_Newton + Dphi     //
// This is the function coupling(a)                        //
//=========================================================//
double coupling_function(double a){
#if defined(FOFRGRAVITY) || defined(MBETAMODEL)

  return 2.0 * beta_of_a(a) * beta_of_a(a);

#elif defined(DGPGRAVITY)

  return 1.0 / (3.0 * beta_DGP(a));

#elif defined(KMOFLAGE)

  // PofX_ofa not implemented...
  return 2.0 * beta * beta / PofX_ofa(a);

#elif defined(BRANSDICKE)

  return JBD_GeffG_of_a(a) - 1.0;

#else

  // Define your function here...
  return 0.0;

#endif
}

//=========================================================//
// The time-dependent effective Newton's constant          //
// This is sometimes called mu(a,k)                        //
// k is assumed to be in units of h/Mpc                    //
//=========================================================//
double GeffoverG(double a, double k){
  double mu = 1.0; 
  
  if(! modified_gravity_active ) return mu;
  if(  use_lcdm_growth_factors ) return mu;

#if defined(FOFRGRAVITY) || defined(MBETAMODEL)

  double mass2a2 = a * a * mass2_of_a(a);
  double k2 = pow2(k * INVERSE_H0_MPCH);
  if(k2 == 0) return mu;
  mu *= 1.0 + 2.0 * beta_of_a(a) * beta_of_a(a) * k2 / ( k2 + mass2a2 ); 
  return mu;

#elif defined(DGPGRAVITY)

  mu *= 1.0 + 1.0/(3.0 * beta_DGP(a));
  return mu;

#elif defined(KMOFLAGE)

  mu *= 1.0 + coupling_function(a);
  return mu;

#elif defined(BRANSDICKE)
  
  mu *= 1.0 + coupling_function(a);
  return mu;

#else

  return mu;

#endif
}

//================================================================
// The massive neutrino modifications enters as en effective 
// mu = (f_cdm * delta_cdm + f_nu * delta_nu) / delta_cdm so we treat it here
// using the transfer-function to compute the term delta_nu/delta_cdm
// The two functions below are (mu = 1 for LCDM):
// D1'' - beta (mu * munu_1LPT) D1 =  beta (mu * munu_1LPT) D1
// D2'' - beta (mu * munu_2LPT) D2 = -beta (mu * munu_2LPT)/2 (D2 - D1 * D1)
//================================================================
double GeffoverG_neutrino_1LPT(double a, double k){
  double munu = 1.0;
#ifdef MASSIVE_NEUTRINOS
  if(OmegaNu > 1e-6)
    munu = (Omega - OmegaNu)/Omega + OmegaNu/Omega * get_nu_transfer_function(k,a) / get_cdm_baryon_transfer_function(k,a);
#endif
  return munu;
}

double GeffoverG_neutrino_2LPT(double a, double k){
  double munu = 1.0;
#ifdef MASSIVE_NEUTRINOS
  if(OmegaNu > 1e-6)
    munu = (Omega - OmegaNu)/Omega;
  //  munu = (Omega - OmegaNu)/Omega + OmegaNu/Omega * get_nu_transfer_function(k,a) / get_cdm_baryon_transfer_function(k,a);
#endif
  return munu;
}

//=========================================================//
// Compute the fifth-force                                 //
// Have to be called right after PtoMesh is run as this    //
// routine assumes we have density(x) in mgarray_two       // 
// and density(k) stored in P3D                            //
//=========================================================//
void ComputeFifthForce(){
  timer_start(_ComputeFifthForce);
  if(! modified_gravity_active){
    if(ThisTask == 0) printf("===> Error: Modified gravity solver is not active so we cannot call this routine!\n");
    timer_stop(_ComputeFifthForce);
    return;
  }

#if defined(FOFRGRAVITY) || defined(MBETAMODEL)

  ComputeFifthForce_PotentialScreening();

#elif defined(DGPGRAVITY)

  ComputeFifthForce_DensityScreening();

#elif defined(KMOFLAGE)

  ComputeFifthForce_GradientScreening();

#elif defined(BRANSDICKE)

  ComputeFifthForce_TimeDepGeffModels();

#else

  // Choose your screening type here...

#endif

  if(ThisTask == 0) printf("\n");
  timer_stop(_ComputeFifthForce);
}

//=========================================================//
// The screening factor for chameleon like models          //
//=========================================================//

double screening_factor_potential(double a, double phinewton){
  if( ! include_screening ) return 1.0;
  if( phinewton >= 0.0) return 1.0;

#if defined(FOFRGRAVITY)

  double phicrit = 1.5 * fofr0 * pow((Omega + 4.0 * (1.0 - Omega)) / (1.0 / (a * a * a) * Omega + 4.0 * (1.0 - Omega)), nfofr + 1.0);
  double screenfac = fabs(phicrit / phinewton);
  if(screenfac > 1.0) screenfac = 1.0;
  return screenfac;

#elif defined(MBETAMODEL)

  // Expression for a general m(a), beta(a) model
  double phicrit = phi_of_a(a);
  double screenfac = fabs(phicrit / phinewton);
  if(screenfac > 1.0) screenfac = 1.0;
  return screenfac;

#else

  // Define your function here...
  return 1.0;

#endif
}

//=========================================================//
// The screening factor depending on density               //
// Vainshtein screening                                    //
//=========================================================//

double screening_factor_density(double a, double density){
  if( ! include_screening ) return 1.0;

#if defined(DGPGRAVITY)

  double fac = 8.0 / 9.0 * Omega * pow( rcH0_DGP / beta_DGP(a) , 2) * (1.0 + density);
  if(fac < 1e-5) return 1.0;
  double screenfac = 2.0*(sqrt(1.0 + fac) - 1.0)/fac;
  return screenfac;

#else

  // Define your function here...
  return 1.0;

#endif
}

//=========================================================//
// The screening factor depending on the gradient          //
// of the Newtonian potential (DPhi2 = |DPhi|^2)           //
// For a P(X) lagrangian with a conformal                  //
// coupling e^{beta phi / Mpl} we have                     //
// P_X^2 X = (2beta^2) Mpl^2  |DPhi|^2 determining X       //
// and then screening_factor is then Min[1, 1/P_X ]      //
// In this case the coupling above is 2beta^2 / P_X(a)     //
//=========================================================//

double screening_factor_gradient(double a, double DPhi2){
  if( ! include_screening ) return 1.0;

#if defined(KMOFLAGE)

  // inverse_phi_function not implemented...
  double X = inverse_phi_function(DPhi2);
  dobule screenfac = 1.0 / PofX(X);
  if(screenfac > 1.0) screenfac = 1.0;
  return screenfac;

#else

  // Define your function here...
  return 1.0;

#endif
}

//=========================================================//
// The time-dependent factor "f" in the 2LPT equation      //
//          D2'' - beff D2 = -beff D1^2 * f                //
//                                                         //
// This is useful for DGP-like theories which have         //
// just a time-depdendent modification here                //
//                                                         //
// This function is only really in play if                 //
// SCALEDEPENDENT is not defined                           //
//=========================================================//
double Factor_2LPT(double a){
  if(! modified_gravity_active || use_lcdm_growth_factors) return 1.0;

#if defined(FOFRGRAVITY) || defined(MBETAMODEL)

  // This is scale-dependent for these models so we compute this elsewhere
  // Only relevant if use_lcdm_growth_factors = 1 for which we return 1.0 above
  return 1.0; 

#elif defined(DGPGRAVITY)

  // For DGP this is the same as the combination used when SCALEDEPDENDENT is set, i.e.
  // 1.0 + 2.0 * second_order_kernel(0, 0, 0, 0, a) / (1.5 * Omega * GeffoverG(a, 0) * a);
  double betadgp = beta_DGP(a);
  return 1.0 - (2.0 * rcH0_DGP * rcH0_DGP * Omega) / (9.0 * a * a * a * pow3( betadgp ) * (1 + 1.0 / (3.0 * betadgp) ) );


#else

  // Define your function here...
  return 1.0;

#endif

}

#ifdef SCALEDEPENDENT

#ifdef FOFRGRAVITY
//=========================================================//
// This is Pi/H0^2 as defined in Bose&Koyama needed        //
// to compute the second order kernel in f(R)              //
// Assuming k in units of h/Mpc                            //
//=========================================================//
double fofr_pi_factor(double k, double a){
  double a3   = a*a*a;
  double fac  = (Omega/a3 + 4.0*(1.0 - Omega));
  double fac0 = (Omega    + 4.0*(1.0 - Omega));
  return pow2(k * INVERSE_H0_MPCH / a) + fac0 * pow( fac / fac0, nfofr + 2.0) / (1.0 + nfofr) / fofr0;
}
#endif

//=========================================================//
// This is the integral kernel gamma_2(k,k1,k2,a)          //
// determining the second order growth-factor              //
//       D2'' - beff D2 = -beff D1(k1)D1(k2) x             //
//        x (1+cos^2 + 2a4H2/beff gamma2)                  //
//=========================================================//
double second_order_kernel(double k, double k1, double k2, double costheta, double a){
  if(! modified_gravity_active || use_lcdm_growth_factors) return 0.0;
  
  double gamma2 = 0.0;

#if defined(FOFRGRAVITY)

  double a3     = a*a*a;
  double fac    = Omega/a3 + 4.0*(1.0 - Omega);
  double fac0   = Omega    + 4.0*(1.0 - Omega);

  gamma2 = - (3.0 * (nfofr+2.0) /(12.0 * (1.0+nfofr) * (1.0+nfofr))) * pow2( k * INVERSE_H0_MPCH / (a * hubble(a)) ) * pow2(Omega/a3) * fac0 * pow(fac/fac0, 2.0*nfofr+3.0) / pow2(fofr0);
  gamma2 /= fofr_pi_factor(k,a) * fofr_pi_factor(k1,a) * fofr_pi_factor(k2,a);

#elif defined(DGPGRAVITY)

  // For completeness, for DGP one should use the now SCALEDEPENDENT version.
  // In terms of Factor_2LPT we have second_order_kernel = 0.5 * (Factor_2LPT(a) - 1) * (1.5 * Omega * GeffoverG(a) * a) * ( 1 - cos^2theta )
  gamma2 = -1.0/6.0/pow3(beta_DGP(a)) * pow2(rcH0_DGP / hubble(a)) * pow2(Omega/(a*a*a)) * (1.0 - costheta*costheta);

#elif defined(MBETAMODEL)

  // The general expression for gamma2 in terms of beta and m
  double m2 = mass2_of_a(a);
  double dm2da = dmass2_of_ada(a);
  double pifac_k  = pow2(k  * INVERSE_H0_MPCH / a) + m2;
  double pifac_k1 = pow2(k1 * INVERSE_H0_MPCH / a) + m2;
  double pifac_k2 = pow2(k2 * INVERSE_H0_MPCH / a) + m2;
  gamma2 = (dm2da/m2) * pow2(beta_of_a(a)) * Omega / (2.0 * a * a) * (m2/pifac_k1) * (m2/pifac_k2) * (1.0 - m2/pifac_k) / pow2(hubble(a));

#else

  // ...

#endif
  
  // Common for all models. If test to avoid dividing by zero
  if(k2 > 1e-3 && k1 > 1e-3){
    gamma2 += 1.5 * Omega / (a * a * a * hubble(a) * hubble(a) ) * 
      ( (GeffoverG(a,k) - GeffoverG(a,k1))*(k1/k2) + (GeffoverG(a,k) - GeffoverG(a,k2))*(k2/k1) ) * costheta/2.0;
  }
  return gamma2;
}

#endif

