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

#if defined(FOFRGRAVITY)

  if(ThisTask == 0){
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

  if(ThisTask == 0){
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

  if(ThisTask == 0){
    printf("============================================\n");
    printf("Running with Modified Gravity, (m(a),b(a))  \n");
    printf("============================================\n");
    printf("  Modified gravity active = %i\n", modified_gravity_active);
    printf("  Include screening       = %i\n", include_screening);
    printf("  Use LCDM growthfactor   = %i\n", use_lcdm_growth_factors);
    printf("  Input P(k) is for LCDM  = %i\n", input_pofk_is_for_lcdm);
    printf("  Sigma8 is for LCDM      = %i\n", input_sigma8_is_for_lcdm);
    printf("\n");
    fflush(stdout);
  }

  // For (m,beta) models we need to compute Phi_crit(a) and store in spline
  compute_phi_of_a();

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

#if defined(FOFRGRAVITY) 

  strcpy(tag[(*nt)], "fofr0");
  addr[(*nt)] = &fofr0;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "nfofr");
  addr[(*nt)] = &nfofr;
  id[(*nt)++] = FLOAT;

#elif defined(MBETAMODEL)

  strcpy(tag[(*nt)], "fofr0");
  addr[(*nt)] = &fofr0;
  id[(*nt)++] = FLOAT;

  strcpy(tag[(*nt)], "nfofr");
  addr[(*nt)] = &nfofr;
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

#else

  // Define your parameters here...

#endif
}

#if defined(BRANSDICKE)

void   read_JBD_background_from_file();
double JBD_Hubble(double a);
double JBD_dHubbleda(double a);
double JBD_GeffOverG(double a);

GSL_Spline JBD_GeffSpline;
GSL_Spline JBD_HSpline;
GSL_Spline JBD_dHdaSpline;

#endif

//=========================================================//
// Hubble function and derivative                          //
// NB: if these are modified then what we call LCDM growth //
// factors will also be modified accordingly.              //  
//=========================================================//
double hubble(double a){

#if defined(BRANSDICKE)

  // Approximation: will add full solution here
  double OmegaCC = 1.0 - Omega;
  double deltaOmega = - wBD/(6.0*(1.0 + wBD)*(1.0 + wBD)) + 1.0/(1.0 + wBD);
  return sqrt( pow(a, -1.0/(1.0 + wBD)) * ( Omega/(a*a*a) + OmegaCC + deltaOmega ) / (1.0 + deltaOmega) );

#else

  return sqrt(Omega / (a*a*a) + 1.0 - Omega);

#endif
}

double dhubbleda(double a){

#if defined(BRANSDICKE)

  // Approximation: will add full solution here
  double H = hubble(a);
  double deltaOmega = -wBD/(6.0*(1.0 + wBD)*(1.0 + wBD)) + 1.0/(1.0 + wBD);
  return  -1.0/(2.0 * H * a * (1.0 + deltaOmega)) * ( 3.0 * Omega/(a*a*a) ) -  H/(2.0 * a * (1.0 + wBD));

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

  // This is beta using f(R) as an example
  return 1.0/sqrt(6.0);

#else

  // E.g. for the symmetron we have
  // if(a > assb_symm){
  //   return beta_symm * sqrt(1.0 - pow3(a/assb_symm));
  // } else {
  //   return 0.0;
  // }

  // Define your function here...
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

  // This is m^2/H0^2 using f(R) as an example
  double a3    = a * a * a;
  double fac   = Omega/a3 + 4.0 * (1.0-Omega);
  double fac0  = Omega    + 4.0 * (1.0-Omega);
  double mass2 = fac0 * pow( fac / fac0, nfofr + 2.0) / ((1.0 + nfofr) * fofr0);
  return mass2; 

#else

  // E.g. for the symmetron we have
  // if(a > assb_symm){
  //   return 0.5 * pow2(2998.0 / range_symm ) * (1.0 - pow3(a/assb_symm));
  // } else {
  //   return 1.0; // Unimporant as beta == 0
  // }

  // Define your function here...
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
// Not used for DGP                                        //
//=========================================================//
#if defined(MBETAMODEL)
double phi_of_a(double a){

  // E.g. for the symmetron we have
  // return 3.0 * Omega / pow3(assb_symm) * pow2(range_symm / 2998.0);

  // Use what we computed in src/new_cosmo.h::compute_phi_of_a()
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

  // Approximation: will add full solution here
  return pow(a, -1.0/(1.0 + wBD)) - 1.0;

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
  if(! modified_gravity_active ) return 1.0;
  if(  use_lcdm_growth_factors ) return 1.0;

#if defined(FOFRGRAVITY) || defined(MBETAMODEL)

  double mass2a2 = a * a * mass2_of_a(a);
  double k2 = pow2(k * INVERSE_H0_MPCH);
  if(k2 == 0) return 1.0;
  return 1.0 + 2.0 * beta_of_a(a) * beta_of_a(a) * k2 / ( k2 + mass2a2 ); 

#elif defined(DGPGRAVITY)

  return 1.0 + 1.0/(3.0 * beta_DGP(a));

#elif defined(KMOFLAGE)

  return 1.0 + coupling_function(a);

#elif defined(BRANSDICKE)
  
  return 1.0 + coupling_function(a);

#else

  return 1.0;

#endif
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
// and then screening_function is then Min[1, 1/P_X ]      //
// In this case the coupling above is 2beta^2 / P_X(a)     //
//=========================================================//

double screening_function_gradient(double a, double DPhi2){
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
  if(! modified_gravity_active) return 1.0;

#if defined(FOFRGRAVITY) || defined(MBETAMODEL)

  // This is scale-dependent for these models so we compute this elsewhere
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
  double fac  = (Omega/a3 + 4.0*(1 - Omega));
  double fac0 = (Omega    + 4.0*(1 - Omega));
  return pow2(k * INVERSE_H0_MPCH / a) + fac*fac*fac / (2.0 * fofr0 * fac0 * fac0);
}
#endif

//=========================================================//
// This is the integral kernel gamma_2(k,k1,k2,a)          //
// determining the second order growth-factor              //
//       D2'' - beff D2 = -beff D1(k1)D1(k2) x             //
//        x (1+cos^2 + 2a4H2/beff gamma2)                  //
//=========================================================//
double second_order_kernel(double k, double k1, double k2, double costheta, double a){
#if defined(FOFRGRAVITY)

  double a3     = a*a*a;
  double fac    = Omega/a3 + 4.0*(1 - Omega);
  double gamma2 = - 9.0/48.0 * pow2( k * INVERSE_H0_MPCH / (a * hubble(a)) ) * pow2(Omega/a3) * pow5(fac)/(pow2(fofr0) * pow4(4 - 3*Omega));
  gamma2 /= fofr_pi_factor(k,a) * fofr_pi_factor(k1,a) * fofr_pi_factor(k2,a);
  return gamma2;

#elif defined(DGPGRAVITY)

  // For completeness, for DGP one should use the now SCALEDEPENDENT version.
  // In terms of Factor_2LPT we have second_order_kernel = 0.5 * (Factor_2LPT(a) - 1) * (1.5 * Omega * GeffoverG(a) * a) * ( 1 - cos^2theta )
  double gamma2 = -1.0/6.0/pow3(beta_DGP(a)) * pow2(rcH0_DGP / hubble(a)) * pow2(Omega/(a*a*a)) * (1.0 - costheta*costheta);
  return gamma2;

#else

  return 0.0;

#endif
}

#endif

