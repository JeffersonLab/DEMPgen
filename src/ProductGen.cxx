/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  This class uses GSL MathMore's MultiRootFinder to find the
  momenta of the produced particles. In this case, a proton,
  and a negative pion.
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <vector>
#include <stdio.h>
#include <iostream>

#include "ProductGen.hxx"
#include "CustomRand.hxx"
#include "Constants.hxx"

#include "TMath.h"
#include "TRandom3.h"

#include "json/json.h"

using namespace std;
using namespace constants;

ProductGen::ProductGen(Particle* inInteraction, Particle* inTarget):
  Interaction(inInteraction), Target(inTarget)
{

  extern Json::Value obj;

  char AngleGenName[100] = "AngleGen";
  double dummy[2] = {0,1};
  double ThetaRange[2] = {obj["prod_pion_thetamin"].asDouble()*TMath::DegToRad(),
                          obj["prod_pion_thetamax"].asDouble()*TMath::DegToRad()};
  double PhiRange[2] = {0, 360*TMath::DegToRad()};
  AngleGen = new CustomRand(AngleGenName, dummy,
                            ThetaRange, PhiRange);
  CoinToss = new TRandom3();

  W_in_val = W_in();

  Proton = new Particle();
  Pion = new Particle();
  Final = new Particle();

  Initial = new Particle();
  *Initial = *Interaction+*Target;

  UnitVect = new TVector3(0,0,1);
  F = new TF1("F",
              "[6]-sqrt([7]**2+x**2)-sqrt([8]**2+([3]-[0]*x)**2+([4]-[1]*x)**2+([5]-[2]*x)**2)",
              0, 12000);

}

void ProductGen::SetInteraction(Particle * inInteraction)
{
  Interaction = inInteraction;
  W_in_val = W_in();
}

void ProductGen::SetTarget(Particle * inTarget)
{
  Target = inTarget;
  W_in_val = W_in();
}

double ProductGen::W_in()
{
  return (*Interaction+*Target).Mag2();
}

double ProductGen::W_out()
{
  return (*Proton+*Pion).Mag2();
}

double ProductGen::Qsq_in()
{
  return - Interaction->Mag2();
}

int ProductGen::Solve()
{
  double theta = AngleGen->Theta();
  double phi = AngleGen->Phi();

  return this->Solve(theta, phi);
}

int ProductGen::Solve(double theta, double phi)
{
  *Initial = *Interaction+*Target;

  W_in_val = W_in();

  if (W_in_val<0){
    //cout << "W < 0 " << endl;
    return 1;
  }

  //cout << proton_mass_mev << endl;

  UnitVect->SetTheta(theta);
  UnitVect->SetPhi(phi);
  UnitVect->SetMag(1);

  double pars[9];
  pars[0] = UnitVect->X();
  pars[1] = UnitVect->Y();
  pars[2] = UnitVect->Z();
  pars[3] = Initial->Px();
  pars[4] = Initial->Py();
  pars[5] = Initial->Pz();
  pars[6] = Initial->E();
  pars[7] = pion_mass_mev;
  pars[8] = proton_mass_mev;

  F->SetParameters(pars);

  double P = F->GetX(0, 0, pars[6], 0.0001, 10000);

  //std::cout << "Zero: " << F->Eval(P) << std::endl;

  Particle * Pion1 = new Particle(pion_mass_mev,
                                  P*pars[0],
                                  P*pars[1],
                                  P*pars[2]);
  *Pion = *Pion1;

  Particle * Proton1 = new Particle();
  *Proton1 = *Initial-*Pion;
  *Proton = *Proton1;

  if (TMath::Abs(F->Eval(P)) > 1){
    delete Pion1;
    delete Proton1;
    return 1;
  }

  if (!SolnCheck()){
    delete Pion1;
    delete Proton1;
    return 1;
  }

  //Check for Second solution:
  double P2 = F->GetX(0, P+100, pars[6], 0.0001, 10000);

  if (TMath::Abs(F->Eval(P2))> 1){
    //No second soln
    delete Pion1;
    delete Proton1;
    return 0;
  }

  //Try second solution
  Particle * Pion2 = new Particle(pion_mass_mev,
                                  P*pars[0],
                                  P*pars[1],
                                  P*pars[2]);
  *Pion = *Pion2;

  Particle * Proton2 = new Particle();
  *Proton2 = *Initial - * Pion;
  *Proton = *Proton2;

  if (SolnCheck()){
    //Toss a coin
    if (CoinToss->Uniform(0,1)>0.5){
      delete Pion1;
      delete Pion2;
      delete Proton1;
      delete Proton2;
      return 0; // Keep second solution
    }
  }

  //Either SolnCheck or coin toss failed
  //Revert to original solution

  *Proton = *Proton1;
  *Pion = *Pion1;

  delete Pion1;
  delete Pion2;
  delete Proton1;
  delete Proton2;

  return 0;
}

bool ProductGen::SolnCheck()
{
  // Double Checking for solution viability
  if (TMath::Abs(proton_mass_mev-Proton->M())>1){
    //cerr << "Mass Missmatch" << endl;
    //cerr << TMath::Abs(proton_mass_mev-Proton->M()) << endl;
    return false;
  }
  if (TMath::Abs(W_in()-W_out())>1){
    //cerr << "W Missmatch" << endl;
    //cerr << TMath::Abs(W_in()-W_out()) << endl;
    return false;
  }
  *Final = *Proton + *Pion;

  if (TMath::Abs(Initial->Px()-Final->Px())>1){
    //cerr << "Px Missmatch" << endl;
    //cerr << TMath::Abs(Initial->Px()-Final->Px()) << endl;
    return false;
  }

  if (TMath::Abs(Initial->Py()-Final->Py())>1){
    //cerr << "Py Missmatch" << endl;
    //cerr << TMath::Abs(Initial->Py()-Final->Py()) << endl;
    return false;
  }

  if (TMath::Abs(Initial->Pz()-Final->Pz())>1){
    //cerr << "Pz Missmatch" << endl;
    //cerr << TMath::Abs(Initial->Pz()-Final->Pz()) << endl;
    return false;
  }

  if (TMath::Abs(Initial->E()-Final->E())>1){
    return false;
  }
  return true;
}

Particle * ProductGen::ProdPion()
{
  return Pion;
}

Particle * ProductGen::ProdProton()
{
  return Proton;
}
