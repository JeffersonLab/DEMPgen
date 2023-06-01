# ifndef legacy_H
# define legacy_H

#include <iostream>
#include <fstream>
#include <string>

#include "eic_pim.h"

#include "TString.h"
#include "TF1.h"
#include "tssa_sig_Para.h"

void Exclusive_Pion_Prodoction(pim);
void Exclusive_Omega_Prodoction(pim);

double fSig_fpi_sigT (double, double); 
double fSig_fpi_sigL (double, double);

#endif

