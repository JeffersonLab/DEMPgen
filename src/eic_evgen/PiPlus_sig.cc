// 21/12/22 - SJDK - The PiPlus cross section routine stripped out of the DEMP_Reaction.cc code in Process routine
// Inputs are -t, W, Qsq and epsilon (in GeV where appropriate)

// Note, there isn't really any reason to split out the parameter determination from this, they could be done in the same file. The PiPlus parameters are a mess anyway.

#include "PiPlus_sig.h"
#include "PiPlus_sig_Param.h"
#include "reaction_routine.h"
#include "TF1.h"

double GetPiPlus_CrossSection(double ft, double fw, double fqsq, double feps ){

  double_t sig_total;

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
 
  // -------------------------------------------------------------------------------------------
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
 
  // -------------------------------------------------------------------------------------------
 
  fSig_VR = fSig_T + feps * fSig_L;

  sig_total = fSig_VR;

  return sig_total;
}
