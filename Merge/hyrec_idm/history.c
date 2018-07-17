/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.c: compute reionization history                                               */
/*                                                                                               */
/*         Version: January 2011                                                                 */
/*         Revision history:                                                                     */
/*            - written November 2010                                                            */
/*            - January 2011: changed various switches (notably for post-Saha expansions)        */
/*                             so that they remain valid for arbitrary cosmologies               */
/*************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"

/*************************************************************************************************
Cosmological parameters Input/Output
*************************************************************************************************/

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param) {

  /* Cosmology */
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter CMB temperature today [Kelvin]: ");
  fscanf(fin, "%lg", &(param->T0));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter baryon density, omega_bh2: ");
  fscanf(fin, "%lg", &(param->obh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter total matter (CDM+baryons) density, omega_mh2: ");
  fscanf(fin, "%lg", &(param->omh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter curvature, omega_kh2: ");
  fscanf(fin, "%lg", &(param->okh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter dark energy density, omega_deh2: ");
  fscanf(fin, "%lg", &(param->odeh2));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter dark energy equation of state parameters, w wa: ");
  fscanf(fin, "%lg %lg", &(param->w0), &(param->wa));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter primordial helium mass fraction, Y: ");
  fscanf(fin, "%lg", &(param->Y));
  if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter effective number of neutrino species, N_nu_eff: ");
  fscanf(fin, "%lg", &(param->Nnueff));

  param->nH0 = 11.223846333047*param->obh2*(1.-param->Y);  /* number density of hudrogen today in m-3 */
  param->fHe = param->Y/(1-param->Y)/3.97153;              /* abundance of helium by number */


  /* Redshift range */
  param->zstart = 8000.;
  param->zend = 0.;
  param->dlna = 8.49e-5;
  param->nz = (long) floor(2+log((1.+param->zstart)/(1.+param->zend))/param->dlna);

  if (fout!=NULL && PROMPT==1) fprintf(fout, "\n");
}

/*************************************************************************************
Hubble expansion parameter in sec^-1
*************************************************************************************/

double rec_HubbleConstant(REC_COSMOPARAMS *param, double z) {

   int j;
   double rho, rho_i[5], ainv;
   double ogh2;

   ainv = 1.+z; /* inverse scale factor */

   /* Matter */
   rho_i[0] = param->omh2 *ainv*ainv*ainv;

   /* Curvature */
   rho_i[1] = param->okh2 *ainv*ainv;

   /* Dark energy */
   rho_i[2] = param->odeh2 * pow(ainv, 3*(1+param->w0)) * exp(3*param->wa* (log(ainv)-1.+1./ainv));

   /* Get radiation density.
    * Coefficient based on 1 AU = 1.49597870691e11 m (JPL SSD) */
   ogh2 = 4.48162687719e-7 * param->T0 * param->T0 * param->T0 * param->T0;
   rho_i[3] = ogh2 *ainv*ainv*ainv*ainv;

   /* Neutrino density -- scaled from photons assuming lower temperature */
   rho_i[4] = 0.227107317660239 * rho_i[3] * param->Nnueff;

   /* Total density, including curvature contribution */
   rho = 0.; for(j=0; j<5; j++) rho += rho_i[j];

   /* Conversion to Hubble rate */
   return( 3.2407792896393e-18 * sqrt(rho) );
}

/*****************************************************************************************
Matter temperature -- 1st order steady state, from Hirata 2008
******************************************************************************************/

double rec_Tmss(double xe, double Tr, double H, double fHe, double nH, double energy_rate) {

  double chi_heat;

  //chi_heat = (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski

  // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
  if (xe < 1.) {
    chi_heat = 0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035));
    if (chi_heat > 1.) chi_heat = 1.;
  }
  else
    chi_heat = 1.;

  return Tr/(1.+H/4.91466895548409e-22/Tr/Tr/Tr/Tr*(1.+xe+fHe)/xe)
    +2./3./kBoltz*chi_heat/nH*energy_rate/(4.91466895548409e-22*pow(Tr,4)*xe);

  /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
}

/******************************************************************************************
Matter temperature evolution derivative
******************************************************************************************/

double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe, double nH, double energy_rate, double dQdlna_baryons) {

  double chi_heat;
  double result;

  //chi_heat = (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski
  // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
  if (xe < 1.) {
    chi_heat = 0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035));
    if (chi_heat > 1.) chi_heat = 1.;
  }
  else
    chi_heat = 1.;
  result = -2.*Tm + 4.91466895548409e-22*Tr*Tr*Tr*Tr*xe/(1.+xe+fHe)*(Tr-Tm)/H
    +2./3./kBoltz*chi_heat/nH*energy_rate/(1.+xe+fHe)/H;
  if(dQdlna_baryons != 0.0){
    // printf("result %e dQdlna_baryons %e\n",result,dQdlna_baryons);
    result += dQdlna_baryons;
  }
  return result;
}
/******************************************************************************************
Dark Matter temperature evolution derivative
******************************************************************************************/

double rec_dTdmeffdlna(double Tdmeff,  double dQdlna_dmeff) {
  // printf("adia %e scattering %e\n",Tdmeff,dQdlna_dmeff);
  return -2.*Tdmeff + dQdlna_dmeff;
}

/**********************************************************************************************
Second order integrator using derivative from previous time steps
Evolves xe only, assumes Tm is given by the steady-state solution
***********************************************************************************************/

void rec_get_xe_next1(REC_COSMOPARAMS *param, double z1, double xe_in, double *xe_out,
                      HRATEEFF *rate_table, int func_select, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
		      double **logfminus_hist, double *logfminus_Ly_hist[],
                      double *z_prev, double *dxedlna_prev, double *z_prev2, double *dxedlna_prev2) {

  double dxedlna, Tr, nH, ainv, H, Tm;

    Tr = param->T0 * (ainv=1.+z1);
    nH = param->nH0 * ainv*ainv*ainv;
    H = rec_HubbleConstant(param, z1);
    Tm = rec_Tmss(xe_in, Tr, H, param->fHe, nH*1e-6, energy_injection_rate(param,z1));

    #if (MODEL == PEEBLES)
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                           rec_HPeebles_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1));
    #elif (MODEL == RECFAST)
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                           rec_HRecFast_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1));
    #elif (MODEL == EMLA2s2p)
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                           rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), rate_table);
    #else
        dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
	  func_select==FUNC_H2G  ? rec_HMLA_2photon_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, rate_table, twog_params,
							    param->zstart, param->dlna, logfminus_hist, logfminus_Ly_hist, iz, z1, energy_injection_rate(param,z1))
	  :rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), rate_table);
    #endif

    *xe_out = xe_in + param->dlna * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));

    *z_prev2       = *z_prev;
    *dxedlna_prev2 = *dxedlna_prev;
    *z_prev        = z1;
    *dxedlna_prev  = dxedlna;

}

/**********************************************************************************************
Second order integrator using derivative from previous time steps
Evolves xe and Tm simultaneously
***********************************************************************************************/

void rec_get_xe_next2(REC_COSMOPARAMS *param, double z1, double xe_in, double Tm_in, double Tdmeff_in,double *xe_out, double *Tm_out,
                      double *Tdmeff_out,
                      HRATEEFF *rate_table, int func_select, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
		                  double **logfminus_hist, double *logfminus_Ly_hist[],
                      double *z_prev, double *dxedlna_prev, double *dTmdlna_prev,double *dTdmeffdlna_prev,
                      double *z_prev2, double *dxedlna_prev2, double *dTmdlna_prev2, double *dTdmeffdlna_prev2) {

    double dxedlna, dTmdlna, dTdmeffdlna, Tr, nH, ainv, H;

    Tr = param->T0 * (ainv=1.+z1);
    nH = param->nH0 * ainv*ainv*ainv;
    H = rec_HubbleConstant(param, z1);
    // printf("before m_in  %e Tdmeff_in  %e\n",Tm_in,Tdmeff_in);

    #if (MODEL == PEEBLES)
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                            rec_HPeebles_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1));
    #elif (MODEL == RECFAST)
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                            rec_HRecFast_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1));
    #elif (MODEL == EMLA2s2p)
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                                            rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), rate_table);
    #else
         dxedlna = func_select==FUNC_HEI  ? rec_helium_dxedt(xe_in, param->nH0, param->T0, param->fHe, H, z1)/H:
                   func_select==FUNC_H2G  ? rec_HMLA_2photon_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, rate_table, twog_params,
                                                                     param->zstart, param->dlna, logfminus_hist, logfminus_Ly_hist, iz, z1,
																	 energy_injection_rate(param,z1)):
	           func_select==FUNC_HMLA ? rec_HMLA_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1), rate_table)
                                           :rec_HPeebles_dxedlna(xe_in, nH*1e-6, H, Tm_in*kBoltz, Tr*kBoltz, energy_injection_rate(param,z1)); /* used for z < 20 only */
    #endif

	 dTmdlna = rec_dTmdlna(xe_in, Tm_in, Tr, H, param->fHe, nH*1e-6, energy_injection_rate(param,z1),calculate_dQdnla_baryons(param,z1,xe_in,Tm_in,Tdmeff_in));
	 dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_in, calculate_dQdnla_dmeff(param,z1,xe_in,Tm_in,Tdmeff_in));

    *xe_out = xe_in + param->dlna * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));
    *Tm_out = Tm_in + param->dlna * (1.25 * dTmdlna - 0.25 * (*dTmdlna_prev2));
    *Tdmeff_out = Tdmeff_in + param->dlna * (1.25 * dTdmeffdlna - 0.25 * (*dTdmeffdlna_prev2));
    // printf("after Tdmeff_out  %e Tdmeff_in %e dTdmeffdlna %e dTdmeffdlna_prev2 %e\n",*Tdmeff_out,Tdmeff_in,dTdmeffdlna,*dTdmeffdlna_prev2);

    *z_prev2       = *z_prev;
    *dxedlna_prev2 = *dxedlna_prev;
    *dTmdlna_prev2 = *dTmdlna_prev;
    *dTdmeffdlna_prev2 = *dTdmeffdlna_prev;
    *z_prev        = z1;
    *dxedlna_prev  = dxedlna;
    *dTmdlna_prev  = dTmdlna;
    *dTdmeffdlna_prev  = dTdmeffdlna;

    // printf("end Tm_in %e Tm_out %e Tdmeff_in %e Tdmeff_Out %e\n",Tm_in,*Tm_out,Tdmeff_in,*Tdmeff_out);

}

/****************************************************************************************************
Builds a recombination history
****************************************************************************************************/

void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       double *xe_output, double *Tm_output, double *Tdmeff_output) {


   long iz;
   double **logfminus_hist;
   double *logfminus_Ly_hist[3];
   double H, z,z_m_1, z_prev, dxedlna_prev, z_prev2, dxedlna_prev2, dTmdlna_prev, dTmdlna_prev2, dTdmeffdlna, dTdmeffdlna_prev, dTdmeffdlna_prev2;
   double Delta_xe;
   double dTdlna;
   /* history of photon occupation numbers */
   logfminus_hist = create_2D_array(NVIRT, param->nz);
   logfminus_Ly_hist[0] = create_1D_array(param->nz);   /* Ly-alpha */
   logfminus_Ly_hist[1] = create_1D_array(param->nz);   /* Ly-beta  */
   logfminus_Ly_hist[2] = create_1D_array(param->nz);   /* Ly-gamma */


   z = param->zstart;
   Tdmeff_output[0] =  param->Tdmeff_initial; // param->Tdmeff_initial calculated within thermodynamics.c assuming interaction has no feedback onto baryons.
   /********* He II + III Saha phase *********/
   Delta_xe = 1.;   /* Delta_xe = xHeIII */
   for(iz=0; iz<param->nz && Delta_xe > 1e-9; iz++) {
      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;

      xe_output[iz] = rec_sahaHeII(param->nH0,param->T0,param->fHe,z, &Delta_xe);
      //Crude integration to be checked. We assume the interaction to have no feedback on baryons at these redshifts.
      Tm_output[iz] = param->T0 * (1.+z);
      if(iz==0){
        Tdmeff_output[iz] = param->Tdmeff_initial;
      }
      else if(iz==1 || iz==2){
        dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_output[iz-1],calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-1],Tm_output[iz-1],Tdmeff_output[iz-1]));
        Tdmeff_output[iz] =Tdmeff_output[iz-1]+ param->dlna * dTdmeffdlna;
      }
      else {
        z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-1)) - 1.;
        z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
        dTdmeffdlna_prev2 = rec_dTdmeffdlna(Tdmeff_output[iz-3],calculate_dQdnla_dmeff(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
        dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_output[iz-1],calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-1],Tm_output[iz-1],Tdmeff_output[iz-1]));

        Tdmeff_output[iz] =Tdmeff_output[iz-1]+ param->dlna * (1.25 * dTdmeffdlna - 0.25 * (dTdmeffdlna_prev2));
      }
   }
   /******* He I + II post-Saha phase *********/
   Delta_xe = 0.;     /* Delta_xe = xe - xe(Saha) */
   // printf("here ok 2 %e iz %d Tm_output[iz] %e  Tdmeff_output[iz] %e \n",z,iz,Tm_output[iz],Tdmeff_output[iz]);

   for (; iz<param->nz && Delta_xe < 5e-4; iz++) {
      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;

      xe_output[iz] = xe_PostSahaHe(param->nH0,param->T0,param->fHe, rec_HubbleConstant(param,z), z, &Delta_xe);
      Tm_output[iz] = param->T0 * (1.+z);
      z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-1)) - 1.;
      z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
      dTdmeffdlna_prev2 = rec_dTdmeffdlna(Tdmeff_output[iz-3],calculate_dQdnla_dmeff(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
      dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_output[iz-1],calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-1],Tm_output[iz-1],Tdmeff_output[iz-1]));

      Tdmeff_output[iz] =Tdmeff_output[iz-1]+ param->dlna * (1.25 * dTdmeffdlna - 0.25 * (dTdmeffdlna_prev2));
   }
   // printf("here ok 3 %e Tm_output[iz] %e  Tdmeff_output[iz] %e \n",z,Tm_output[iz-1],Tdmeff_output[iz-1]);

   /****** Segment where we follow the helium recombination evolution, Tm fixed to steady state *******/

   z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
   dxedlna_prev2 = (xe_output[iz-2] - xe_output[iz-4])/2./param->dlna;

   z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-2)) - 1.;
   dxedlna_prev = (xe_output[iz-1] - xe_output[iz-3])/2./param->dlna;

   Delta_xe = 1.;  /* Difference between xe and H-Saha value */

   for(; iz<param->nz && (Delta_xe > 1e-4 || z > 1650.); iz++) {

      rec_get_xe_next1(param, z, xe_output[iz-1], xe_output+iz, rate_table, FUNC_HEI, iz-1, twog_params,
		     logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &z_prev2, &dxedlna_prev2);

      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;

      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0*(1.+z), rec_HubbleConstant(param, z), param->fHe, param->nH0*cube(1.+z), energy_injection_rate(param,z));
      z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-1)) - 1.;
      z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
      dTdmeffdlna_prev2 = rec_dTdmeffdlna(Tdmeff_output[iz-3],calculate_dQdnla_dmeff(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
      dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_output[iz-1],calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-1],Tm_output[iz-1],Tdmeff_output[iz-1]));

      Tdmeff_output[iz] =Tdmeff_output[iz-1]+ param->dlna * (1.25 * dTdmeffdlna - 0.25 * (dTdmeffdlna_prev2));
      /* Starting to populate the photon occupation number with thermal values */
      update_fminus_Saha(logfminus_hist, logfminus_Ly_hist, xe_output[iz], param->T0*(1.+z)*kBoltz,
                         param->nH0*cube(1.+z)*1e-6, twog_params, param->zstart, param->dlna, iz, z, 0);

      Delta_xe = fabs(xe_output[iz]- rec_saha_xe_H(param->nH0, param->T0, z));
    }
    // printf("here ok 4 %e Tm_output[iz] %e  Tdmeff_output[iz] %e \n",z,Tm_output[iz-1],Tdmeff_output[iz-1]);


   /******* Hydrogen post-Saha equilibrium phase *********/
   Delta_xe = 0.;  /*Difference between xe and Saha value */

   for(; iz<param->nz && Delta_xe < 5e-5; iz++) {
      z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;

      H = rec_HubbleConstant(param,z);
      xe_output[iz] =  xe_PostSahaH(param->nH0*cube(1.+z)*1e-6, H, kBoltz*param->T0*(1.+z), rate_table, twog_params,
				    param->zstart, param->dlna, logfminus_hist, logfminus_Ly_hist, iz, z, &Delta_xe, MODEL, energy_injection_rate(param,z));
      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0*(1.+z), rec_HubbleConstant(param, z), param->fHe, param->nH0*cube(1.+z), energy_injection_rate(param,z));
      z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-1)) - 1.;
      z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
      dTdmeffdlna_prev2 = rec_dTdmeffdlna(Tdmeff_output[iz-3],calculate_dQdnla_dmeff(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
      dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_output[iz-1],calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-1],Tm_output[iz-1],Tdmeff_output[iz-1]));

      Tdmeff_output[iz] =Tdmeff_output[iz-1]+ param->dlna * (1.25 * dTdmeffdlna - 0.25 * (dTdmeffdlna_prev2));
    }
    // printf("here ok 5 %e Tm_output[iz] %e  Tdmeff_output[iz] %e \n",z,Tm_output[iz-1],Tdmeff_output[iz-1]);

    /******* Segment where we follow the hydrogen recombination evolution with two-photon processes
             Tm fixed to steady state ******/

    z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
    dxedlna_prev2 = (xe_output[iz-2] - xe_output[iz-4])/2./param->dlna;

    z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-2)) - 1.;
    dxedlna_prev = (xe_output[iz-1] - xe_output[iz-3])/2./param->dlna;

    for(; iz<param->nz && 1.-Tm_output[iz-1]/param->T0/(1.+z) < 5e-4 && z > 950.; iz++) { //In original hyrec, z > 700; if scattering with DM is switch on we want to make sure that the full evolution is followed early enough

       rec_get_xe_next1(param, z, xe_output[iz-1], xe_output+iz, rate_table, FUNC_H2G, iz-1, twog_params,
               	      logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &z_prev2, &dxedlna_prev2);
       z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
       // printf("z %e\n",z);
       Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0*(1.+z), rec_HubbleConstant(param, z), param->fHe, param->nH0*cube(1.+z)*1e-6, energy_injection_rate(param,z));
       z_prev = (1.+param->zstart)*exp(-param->dlna*(iz-1)) - 1.;
       z_prev2 = (1.+param->zstart)*exp(-param->dlna*(iz-3)) - 1.;
       dTdmeffdlna_prev2 = rec_dTdmeffdlna(Tdmeff_output[iz-3],calculate_dQdnla_dmeff(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
       dTdmeffdlna = rec_dTdmeffdlna(Tdmeff_output[iz-1],calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-1],Tm_output[iz-1],Tdmeff_output[iz-1]));

       Tdmeff_output[iz] =Tdmeff_output[iz-1]+ param->dlna * (1.25 * dTdmeffdlna - 0.25 * (dTdmeffdlna_prev2));
     }
    // printf("here ok 6 %e Tm_output[iz] %e  Tdmeff_output[iz] %e \n",z,Tm_output[iz-1],Tdmeff_output[iz-1]);

   /******* Segment where we follow the hydrogen recombination evolution with two-photon processes
            AND Tm evolution ******/

    dTmdlna_prev2 = rec_dTmdlna(xe_output[iz-3], Tm_output[iz-3], param->T0*(1.+z_prev2),
                                rec_HubbleConstant(param, z_prev2), param->fHe, param->nH0*cube(1.+z_prev2), energy_injection_rate(param,z_prev2),calculate_dQdnla_baryons(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
    dTmdlna_prev  = rec_dTmdlna(xe_output[iz-2], Tm_output[iz-2], param->T0*(1.+z_prev),
                                rec_HubbleConstant(param, z_prev), param->fHe, param->nH0*cube(1.+z_prev), energy_injection_rate(param,z_prev),calculate_dQdnla_baryons(param,z_prev,xe_output[iz-2],Tm_output[iz-2],Tdmeff_output[iz-2]));
    dTdmeffdlna_prev2 = rec_dTdmeffdlna(Tdmeff_output[iz-3],calculate_dQdnla_dmeff(param,z_prev2,xe_output[iz-3],Tm_output[iz-3],Tdmeff_output[iz-3]));
    dTdmeffdlna_prev = rec_dTdmeffdlna(Tdmeff_output[iz-2], calculate_dQdnla_dmeff(param,z_prev,xe_output[iz-2],Tm_output[iz-2],Tdmeff_output[iz-2]));

    // printf("ready for the loop %e %e %e %e\n",dTmdlna_prev,dTmdlna_prev2,dTdmeffdlna_prev,dTdmeffdlna_prev2);
    for(; iz<param->nz && z > 700.; iz++) {
        rec_get_xe_next2(param, z, xe_output[iz-1], Tm_output[iz-1],Tdmeff_output[iz-1], xe_output+iz, Tm_output+iz,Tdmeff_output+iz, rate_table, FUNC_H2G,
                        iz-1, twog_params, logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &dTmdlna_prev, &dTdmeffdlna_prev,
                        &z_prev2, &dxedlna_prev2, &dTmdlna_prev2, &dTdmeffdlna_prev2);
        z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
     }
     // printf("here ok 7 %e\n",z );

    /***** Segment where we follow Tm as well as xe *****/
    /* Radiative transfer effects switched off here */
    for(; iz<param->nz && z>20.; iz++) {

        rec_get_xe_next2(param, z, xe_output[iz-1], Tm_output[iz-1],Tdmeff_output[iz-1], xe_output+iz, Tm_output+iz,Tdmeff_output+iz, rate_table, FUNC_HMLA,
                        iz-1, twog_params, logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &dTmdlna_prev, &dTdmeffdlna_prev,
                        &z_prev2, &dxedlna_prev2, &dTmdlna_prev2, &dTdmeffdlna_prev2);
        z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
    }
    // printf("here ok 8 %e\n",z );

    /*** For z < 20 use Peeble's model. The precise model does not metter much here as
            1) the free electron fraction is basically zero (~1e-4) in any case and
            2) the universe is going to be reionized around that epoch
         Tm is still evolved explicitly ***/
    for(; iz<param->nz; iz++) {

        rec_get_xe_next2(param, z, xe_output[iz-1], Tm_output[iz-1],Tdmeff_output[iz-1], xe_output+iz,Tm_output+iz, Tdmeff_output+iz, rate_table, FUNC_PEEBLES,
                        iz-1, twog_params, logfminus_hist, logfminus_Ly_hist, &z_prev, &dxedlna_prev, &dTmdlna_prev, &dTdmeffdlna_prev,
                        &z_prev2, &dxedlna_prev2, &dTmdlna_prev2, &dTdmeffdlna_prev2);
        z = (1.+param->zstart)*exp(-param->dlna*iz) - 1.;
    }
    // printf("success in hyrec\n");
    /* Cleanup */
    free_2D_array(logfminus_hist, NVIRT);
    free(logfminus_Ly_hist[0]);
    free(logfminus_Ly_hist[1]);
    free(logfminus_Ly_hist[2]);

}

double onthespot_injection_rate(REC_COSMOPARAMS *param,
			     double z) {

  double annihilation_at_z;
  double rho_cdm_today;
  double u_min;
  double erfc;

  /*redshift-dependent annihilation parameter*/

  if (z>param->annihilation_zmax) {

    annihilation_at_z = param->annihilation*
      exp(-param->annihilation_variation*pow(log((param->annihilation_z+1.)/(param->annihilation_zmax+1.)),2));
  }
  else if (z>param->annihilation_zmin) {

    annihilation_at_z = param->annihilation*
      exp(param->annihilation_variation*(-pow(log((param->annihilation_z+1.)/(param->annihilation_zmax+1.)),2)
					 +pow(log((z+1.)/(param->annihilation_zmax+1.)),2)));
  }
  else {

    annihilation_at_z = param->annihilation*
      exp(param->annihilation_variation*(-pow(log((param->annihilation_z+1.)/(param->annihilation_zmax+1.)),2)
					 +pow(log((param->annihilation_zmin+1.)/(param->annihilation_zmax+1.)),2)));
  }

  rho_cdm_today = param->omh2*1.44729366e-9; /* energy density in Kg/m^3 */

  u_min = (1+z)/(1+param->annihilation_z_halo);

  erfc = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4);

  return (pow(rho_cdm_today,2)/2.99792458e8/2.99792458e8*pow((1.+z),3)*
    (pow((1.+z),3)*annihilation_at_z+param->annihilation_f_halo*erfc)
    +rho_cdm_today*pow((1+z),3)*param->decay)/1.e6/1.60217653e-19;
  /* energy density rate in eV/cm^3/s (remember that annihilation_at_z is in m^3/s/Kg and decay in s^-1) */
  /* note that the injection rate used by recfast, defined in therodynamics.c, is in J/m^3/s. Here we multiplied by 1/1.e6/1.60217653e-19 to convert to eV and cm. */

}

double energy_injection_rate(REC_COSMOPARAMS *param,
					double z) {

  double zp,dz;
  double integrand,first_integrand;
  double factor,result;

  if (param->annihilation > 0.) {

    if (param->has_on_the_spot == 0) {

      /* factor = c sigma_T n_H(0) / H(0) (dimensionless) */
      factor = 2.99792458e8 * 6.6524616e-29 * param->nH0 / (3.2407792896393e-18 * sqrt(param->omh2));

      /* integral over z'(=zp) with step dz */
      dz=1.;

      /* first point in trapezoidal integral */
      zp = z;
      first_integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot_injection_rate(param,zp); // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 8 and 7.5
      result = 0.5*dz*first_integrand;

      /* other points in trapezoidal integral */
      do {

	zp += dz;
	integrand = factor*pow(1+z,8)/pow(1+zp,7.5)*exp(2./3.*factor*(pow(1+z,1.5)-pow(1+zp,1.5)))*onthespot_injection_rate(param,zp); // beware: versions before 2.4.3, there were rwrong exponents: 6 and 5.5 instead of 8 and 7.5
	result += dz*integrand;
	//moment += dz*integrand*(zp-z);

      } while (integrand/first_integrand > 0.02);

      /* test lines for printing energy rate rescaled by (1=z)^6 in J/m^3/s w/o approximation */
      /*                                                                                                                                                                                                                    fprintf(stdout,"%e  %e  %e\n",                                                                                                                                                                                     1.+z,                                                                                                                                                                                                              result/pow(1.+z,6)*1.602176487e-19*1.e6,                                                                                                                                                                           onthespot_injection_rate(param,z)/pow(1.+z,6)*1.602176487e-19*1.e6);                                                                                                                                            */

    }
    else {
      result = onthespot_injection_rate(param,z);
    }

    /* effective energy density rate in eV/cm^3/s  */
    return result;
  }
  else {
    return 0.;
  }

}


double calculate_dQdnla_dmeff(REC_COSMOPARAMS *param,double z,double xe, double Tm,double Tdmeff){
  double a = 1.0/(1+z);
  double dQdlnadmeff;
  double dkappaT_dmeff = calculate_dmeff_rate(param,z,xe,Tm,Tdmeff);
  double H = rec_HubbleConstant(param, z);//in s^-1
  // printf("a %e dkappaT_dmeff %e H %e Tm %e Tdmeff %e\n",a,dkappaT_dmeff,H,Tm,Tdmeff);
  dQdlnadmeff = dkappaT_dmeff * (Tm-Tdmeff)/H/a*2; //Compared to recfast, there is an additional factor 1/a since we calculate dQdlna=adQda=a*(dQdz*dzda) = -dQdz/a
  //Found a factor 2 difference...
  return dQdlnadmeff;
}

double calculate_dQdnla_baryons(REC_COSMOPARAMS *param,double z,double xe, double Tm,double Tdmeff){
  double a = 1.0/(1+z);
  double dQdlna_baryons;
  double dkappaT_dmeff = calculate_dmeff_rate(param,z,xe,Tm,Tdmeff);
  double H = rec_HubbleConstant(param, z);//in s^-1
  double mu = mH/(1. + (1./_not4_ - 1.) * param->Y + xe * (1.-param->Y))* _eV_ / (_c_ * _c_); //mH originally in eV/c2, we convert to Kg
  double rho_dmeff = param->omdmeffh2*cube(1+z); //units are unimportant we take the ratio
  double rho_baryon=param->obh2*cube(1+z);
  // printf("(mu/param->m_dmeff) %e1. + (1./_not4_ - 1.) * param->Y + xe * (1.-param->Y) %e  Y %e xe %e mu %e \n", (mu/param->m_dmeff),1. + (1./_not4_ - 1.) * param->Y + xe * (1.-param->Y), param->Y,xe,mu);
  if(param->dmeff_backreaction_on_baryons == 1){
    dQdlna_baryons = dkappaT_dmeff * (Tdmeff-Tm)/H/a; //Compared to recfast, there is an additional factor 1/a since we calculate dQdlna=adQda=a*(dQdz*dzda) = -dQdz/a
    dQdlna_baryons *= 2.*(mu/param->m_dmeff)*(rho_dmeff/rho_baryon); //From energy conservation
  }
  else{
    dQdlna_baryons = 0.0;
  }
  return dQdlna_baryons;
}

double calculate_dmeff_rate(REC_COSMOPARAMS *param,double z,double xe, double Tm,double Tdmeff){

  double  a, H, rho_baryon, dTdtau_T;
  double mH_kg = mH*(_eV_ / (_c_ * _c_));//convert from eV to kg
  double rate, cn, avg_term,vrms2,result;
  // printf("sigma %e\n", param->sigma_dmeff);
  /* Check if there is any interaction at all. */
  if(param->sigma_dmeff == 0.0){
    return  0.0;
  }

  /* Extract background quantities. */
  a = 1/(1+z);
  H = rec_HubbleConstant(param, z);//in s^-1
  /* baryon density (which includes He), converted into units of kg / m^2 / Mpc */
  // pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  // pba->H0 = param1 * 1.e3 / _c_;
  // pba->h = param1 / 100.;
  rho_baryon = param->obh2* (1.e5 / _c_)* (1.e5 / _c_)  *cube(1+z) * 3.*_c_*_c_/(8.*3.14157*_G_*_Mpc_over_m_);
  // rho_baryon = param->obh2*cube(1+z)*1.44729366e-9*_Mpc_over_m_; /* energy density in Kg/m^2/Mpc */
  // if(z<1000 && z>500)printf("param->obh2* square(1.e5 / _c_) *cube(1+z) %e rho_baryon %e z %e 3.*_c_*_c_/(8.*3.14157*_G_) %e \n",param->obh2* square(1.e5 / _c_) *cube(1+z),rho_baryon,z,3.*_c_*_c_/(8.*3.14157*_G_));
  cn = pow(2.,(param->npow_dmeff + 5.)/2.) * tgamma(3. + param->npow_dmeff / 2.) / (3.*_SQRT_PI_);

  rate = a * (1.-param->Y) * rho_baryon * cn * param->sigma_dmeff/(param->m_dmeff + mH_kg)*(_c_/_Mpc_over_m_); //need to convert from 1/Mpc to 1/s so we multiply by _Mpc_over_m_*_c_

  /* Thermal contribution. */
    if(z>1000){
      vrms2 = param->Vrel_dmeff_rate*param->Vrel_dmeff_rate; //(30e3m/s)^2, see e.g. arXiv:1005.2416
    }
    else{
      vrms2 = param->Vrel_dmeff_rate*param->Vrel_dmeff_rate*pow((1+z)/1000,2);
    }

  // printf("z %e Tb/_m_H_  %e vrms2 %e  param->npow_dmeff %e\n",z, Tm/mH_kg*_k_B_/(_c_*_c_),vrms2/3/(_c_*_c_),param->npow_dmeff);
  avg_term = ((Tm/mH_kg + Tdmeff/param->m_dmeff)*_k_B_+vrms2/3)/(_c_*_c_);
  rate *= pow(avg_term, (param->npow_dmeff + 1.)/2.);

  /*VG: turn off rate artificially at z<100*/
  if(z<param->z_rate_off_dmeff){
    rate = 0;
    }

  /* Save rates */
  result = rate * param->m_dmeff / (param->m_dmeff + mH_kg);

  /* Also fill in c^2 for dmeff will be set later. */

  return result;
}
