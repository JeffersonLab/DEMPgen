// 31/01/23 - SJDK
// New file to get the KPlus cross section via the scaling method that Ali Usman utilised in preliminary kaon studies
// This method relies upon getting the pion cross section before subsequently scaling it

#include "KPlus_sig_Scaling.h"
#include "PiPlus_sig.h"
#include "PiPlus_sig_Param.h"
#include "reaction_routine.h"
#include "TF1.h"
#include "TString.h"

double GetKPlus_CrossSection_Scaling(double ft, double fw, double fqsq, double feps, double meson_mass, TString recoil_hadron){

  double sig_total;

  // --------------------------------------------------------------------------------------------------
  // CKY sigma L and T starts
  // --------------------------------------------------------------------------------------------------
  double lpar0 = 0., lpar1 = 0., lpar2 = 0., lpar3 = 0., lpar4 = 0., lpar5 = 0., lpar6 = 0.;
  double tpar0 = 0., tpar1 = 0., tpar2 = 0., tpar3 = 0., tpar4 = 0.;

  lpar0 = 0.;    lpar1 = 0.;    lpar2 = 0.;    lpar3 = 0.;    lpar4 = 0.;    lpar5 = 0.;    lpar6 = 0.;
  tpar0 = 0.;    tpar1 = 0.;    tpar2 = 0.;    tpar3 = 0.;    tpar4 = 0.;
 
  fSig_L = 0;
  fSig_T = 0;
 
  if ( ( ft > 0. ) && ( ft < 0.15 ) ) {
    PiPlus_sigmaL_Param( fw,  fqsq, lpar0, lpar1, lpar2 , lpar3 , lpar4 , lpar5 , lpar6 );
    TF1 *fitCKYLonglandau = new TF1("sigmaL","landau", 0.0 , 0.15 );
    fitCKYLonglandau->FixParameter( 0 , lpar0 );
    fitCKYLonglandau->FixParameter( 1 , lpar1 );
    fitCKYLonglandau->FixParameter( 2 , lpar2 );
    fSig_L = fitCKYLonglandau->Eval(ft);
    if ( lpar0 == 0 || lpar1 == 0 || lpar2 == 0 )
      fSig_L = 0;
    fitCKYLonglandau = NULL;
    delete fitCKYLonglandau;
  }
  else if ( ( ft > 0.15 ) && ( ft < 0.5 ) ) {
    PiPlus_sigmaL_Param( fw,  fqsq, lpar0, lpar1, lpar2 , lpar3 , lpar4 , lpar5 , lpar6 );
    TF1 *fitCKYLongexpo1 = new TF1("sigmaL","expo", 0.15 , 0.5 );
    fitCKYLongexpo1->FixParameter( 0 , lpar3 );
    fitCKYLongexpo1->FixParameter( 1 , lpar4 );
    fSig_L = fitCKYLongexpo1->Eval(ft);
    if ( lpar3 == 0 || lpar4 == 0 )
      fSig_L = 0;
    fitCKYLongexpo1 = NULL;
    delete fitCKYLongexpo1;
  }
  else if ( ( ft > 0.5 ) && ( ft < 1.3 ) ) {
    PiPlus_sigmaL_Param( fw,  fqsq, lpar0, lpar1, lpar2 , lpar3 , lpar4 , lpar5 , lpar6 );
    TF1 *fitCKYLongexpo2 = new TF1("sigmaL","expo", 0.5 , 1.3 );
    fitCKYLongexpo2->FixParameter( 0 , lpar5 );
    fitCKYLongexpo2->FixParameter( 1 , lpar6 );
    fSig_L = fitCKYLongexpo2->Eval(ft);
    if ( lpar5 == 0 || lpar6 == 0 )
      fSig_L = 0;
    fitCKYLongexpo2 = NULL;
    delete fitCKYLongexpo2;
  }
  else {
    fSig_L = 0;
  }
   // SJDK - 02/06/22 - The validity range here was inconsistent, this only went from 0.0 to 0.15, leaving a gap between 0.15 to 0.2
  // I changed the range to remove this gap. 
  if ( ( ft > 0.0 ) && ( ft < 0.2 ) ) {
    PiPlus_sigmaT_Param( fw,  fqsq, tpar0, tpar1, tpar2 , tpar3 , tpar4 );
    TF1 *fitCKYTranspol2 = new TF1("sigmaL","pol2", 0.0 , 0.2 );
    fitCKYTranspol2->FixParameter( 0 , tpar0 );
    fitCKYTranspol2->FixParameter( 1 , tpar1 );
    fitCKYTranspol2->FixParameter( 2 , tpar2 );
    fSig_T = fitCKYTranspol2->Eval(ft);
    if ( tpar0 == 0 || tpar1 == 0 || tpar2 == 0 )
      fSig_T = 0;
    fitCKYTranspol2 = NULL;
    delete fitCKYTranspol2;
  }
  else if ( ( ft > 0.2 ) && ( ft < 1.3 ) ) {
    PiPlus_sigmaT_Param( fw,  fqsq, tpar0, tpar1, tpar2 , tpar3 , tpar4 );
    TF1 *fitCKYTransexpo = new TF1("sigmaL","expo", 0.2 , 1.3 );
    fitCKYTransexpo->FixParameter( 0 , tpar3 );
    fitCKYTransexpo->FixParameter( 1 , tpar4 );
    fSig_T = fitCKYTransexpo->Eval(ft);
    if ( tpar3 == 0 || tpar4 == 0 )
      fSig_T = 0;
    fitCKYTransexpo = NULL;
    delete fitCKYTransexpo;
  }
 
  // ------------------------------------------------------------------------------------------------
  // Improving the Kaon Sigma_L following GH's fortran code (pole_ration.f) located in main driectory
  // Added by AU on July 10, 2021
  // ------------------------------------------------------------------------------------------------

  double mkg = 0, mpig = 0, hbarc = 0, gpoleKL = 0, gpoleKS = 0, gpolepi = 0, gpoleKhyp = 0, r2_dip = 0, r2_mono = 0, Fpi_mono = 0, Fpi_dip = 0, pmono = 0, Fpi_fit = 0, fpisq = 0, Fkk = 0, fksq = 0, lNpi = 0, lNk = 0, gkhypn = 0, gpinn = 0, dl_poleKhyp = 0, dl_polepi = 0, ratio_Khyp_pi = 0, FT_GeV_neg = 0;

  gpoleKL = -13.3;
  gpoleKS = -3.5; 
  gpolepi = 13.1;

  r2_dip = 0.411;
  r2_mono = 0.431;                                                                                                                                                                                                              

  hbarc = 0.197;
  mkg = meson_mass;
  mpig = 0.13957;

  FT_GeV_neg = -1.0 * ft;

  if ( recoil_hadron == "Lambda" ) {
    gpoleKhyp = gpoleKL;
  }
  else if ( recoil_hadron == "Sigma0" ) {
    gpoleKhyp = gpoleKS;
  }

  Fpi_mono = 1.0 / ( 1.0 + ( r2_mono * fqsq) / (6 * (hbarc * hbarc)));   
  Fpi_dip = 1.0 / ( 1.0 + ( r2_dip * fqsq) / (12 * (hbarc * hbarc))) * ( 1.0 + ( r2_dip * fqsq) / (12 * (hbarc * hbarc)));
  pmono = 0.85;
  
  Fpi_fit = pmono * Fpi_mono + (1 - pmono) * Fpi_dip;
  fpisq = Fpi_fit * Fpi_fit;

  Fkk = 0.9 / (1.0 + fqsq / 0.462 );
  fksq = Fkk * Fkk;

  lNpi = 0.44;
  lNk = (0.44+0.80)/2;

  gkhypn = gpoleKhyp * ((lNk * lNk) - (mkg * mkg)) / ((lNk * lNk) - FT_GeV_neg); 
  gpinn = gpolepi * ((lNpi * lNpi) - (mpig * mpig)) / ((lNpi * lNpi) - FT_GeV_neg);

  dl_poleKhyp = ((gkhypn * gkhypn) * fksq) / ((FT_GeV_neg - (mkg * mkg))*(FT_GeV_neg - (mkg * mkg)));  
  dl_polepi = ((gpinn * gpinn) * fpisq) / ((FT_GeV_neg - (mpig * mpig))*(FT_GeV_neg - (mpig * mpig)));
  
  ratio_Khyp_pi = dl_poleKhyp / dl_polepi; 

  // --------------------------------------------------------------------------------

  fSig_VR = (0.1* fSig_T) + feps * (ratio_Khyp_pi* fSig_L);
  sig_total = fSig_VR;

  return sig_total;
}
