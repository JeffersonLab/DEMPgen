#include "CustomRand.hxx"

#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace std;

CustomRand::CustomRand(const char * fname,
                       double ERange[2],
                       double ThetaRange[2],
                       double PhiRange[2])
{
  extern int Gen_seed;
	
  PartSet = true;

  int seed = 0;
  
  //cout << "Seed setting: " << Gen_seed << endl;
  
//  exit(0);

  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  gRandom -> SetSeed(seed);

  char PhiName[100] = "Phi_";
  strcat(PhiName, fname);
  PhiRand = new TF1(PhiName, "1", PhiRange[0], PhiRange[1]);

  char ThetaName[100] = "Theta_";
  strcat(ThetaName, fname);

  // ThetaRand = new TF1(ThetaName, "1",
  //                     TMath::Cos(ThetaRange[0]),
  //                     TMath::Cos(ThetaRange[1]));

  ThetaRand = new TF1(ThetaName, "0.5*sin(x)", ThetaRange[0], ThetaRange[1]);

  char EName[100] = "E_%c";
  std::sprintf(EName, EName, fname);
  ERand = new TF1(EName, "1", ERange[0], ERange[1]);
}


/*--------------------------------------------------*/

CustomRand::CustomRand(const char * fname,
                       double ERange[2],
                       double ThetaRange[2],
                       double PhiRange[2], 
					   int demp_seed)
{
  PartSet = true;

  int seed = demp_seed;

  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  gRandom -> SetSeed(seed);

  char PhiName[100] = "Phi_";
  strcat(PhiName, fname);
  PhiRand = new TF1(PhiName, "1", PhiRange[0], PhiRange[1]);

  char ThetaName[100] = "Theta_";
  strcat(ThetaName, fname);

  // ThetaRand = new TF1(ThetaName, "1",
  //                     TMath::Cos(ThetaRange[0]),
  //                     TMath::Cos(ThetaRange[1]));

  ThetaRand = new TF1(ThetaName, "0.5*sin(x)", ThetaRange[0], ThetaRange[1]);

  char EName[100] = "E_%c";
  std::sprintf(EName, EName, fname);
  ERand = new TF1(EName, "1", ERange[0], ERange[1]);
}

/*--------------------------------------------------*/

double CustomRand::Phi()
{
  return PhiRand->GetRandom();
}
double CustomRand::Theta()
{
  return ThetaRand->GetRandom();

}
double CustomRand::E()
{
  return ERand->GetRandom();
}
