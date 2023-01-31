# ifndef EIC_H
# define EIC_H

#include <iostream>
#include "TRandom.h"
#include "TRandom2.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>

#include "eic_pim.h"

#include "PiPlus_sig.h"
#include "PiPlus_sig_Param.h"
//#include "KPlusLambda_sig_Param.h"
//#include "KPlusSigma_sig_Param.h"
//#include "Pi0_sig_Param.h"

#include "reaction_routine.h"

#include "json/json.h"
#include "json/json-forwards.h"

using std::vector;

void eic();
//void eic(int, int, int, TString, int, TString);

// 18/01/23 - SJDK- This function is never used since eic() is only called with a json object as the argument. Commented out for now, delete later?
//void eic(int, int, int, TString, int, TString, TString, TString, TString, double, double);
void eic(Json::Value);

extern int fSeed;

void SetEICSeed(int);

TString ExtractParticle(TString);
TString ExtractCharge(TString);

vector<vector<vector<vector<double>>>> ReadCrossSectionPar(TString particle, TString hadron);

#endif



