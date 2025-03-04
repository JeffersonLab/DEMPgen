#include <vector>
#include "json/json.h"

#include "TMath.h"

#include "DEMPEvent.hxx"
#include "Asymmetry.hxx"

#include "SigmaL.hxx"
#include "SigmaT.hxx"
#include "SigmaLT.hxx"
#include "SigmaTT.hxx"
#include "Constants.hxx"

#include "correction.h"

#include "SigmaCalc.hxx"

using namespace std;
using namespace constants;
using namespace TMath;

extern Json::Value obj;


SigmaCalc::SigmaCalc(DEMPEvent* in_VertEvent,
                     DEMPEvent* in_CofMEvent,
                     DEMPEvent* in_RestEvent,
                     DEMPEvent* in_TConEvent):
  VertEvent(in_VertEvent),
  CofMEvent(in_CofMEvent),
  RestEvent(in_RestEvent),
  TConEvent(in_TConEvent)
{

  vector<double> qsq;
  qsq.push_back(4.107);
  qsq.push_back(4.335);
  qsq.push_back(4.845);
  qsq.push_back(5.138);
  qsq.push_back(5.557);
  qsq.push_back(5.948);
  qsq.push_back(6.049);
  qsq.push_back(6.393);
  qsq.push_back(6.778);
  qsq.push_back(6.894);
  qsq.push_back(7.617);
  
  Asyms = new vector<Asymmetry*>(5);
  
  Asyms->at(0) =
    new Asymmetry("asy",
                  "[0]*exp([1]*x)+(-[2]-[0])*exp([3]*x)+[2]",
                  qsq, false);

  Asyms->at(1) =
    new Asymmetry("asysfi",
                  "[0]*exp([1]*x)+[2]",
                  qsq, false);

  Asyms->at(2) =
    new Asymmetry("asy2fi",
                  "[0]*exp([1]*x)+[2]",
                  qsq, false);

  Asyms->at(3) =
    new Asymmetry("asyfpfs",
                  "[0]*exp([1]*x)+[2]",
                  qsq, false);

  Asyms->at(4) =
    new Asymmetry("asy3f",
                  "[0]*exp([1]*x)+[2]",
                  qsq, false);

}

double SigmaCalc::sigma_l()
{
  if(*VertEvent->t_GeV < -1.5)
    return 1e-6;
  else{
    return MySigmaL(*VertEvent->qsq_GeV,
                    -*VertEvent->t_GeV,
                    *VertEvent->w_GeV);
  }
  return MySigmaL(*VertEvent->qsq_GeV,
                  -*VertEvent->t_GeV,
                  *VertEvent->w_GeV);
}

double SigmaCalc::sigma_t()
{
  if (*VertEvent->t_GeV < -1.5)
    return 1e-6;
  else{
    return MySigmaT(*VertEvent->qsq_GeV,
                    -*VertEvent->t_GeV,
                    *VertEvent->w_GeV);
  }
  return MySigmaT(*VertEvent->qsq_GeV,
                  -*VertEvent->t_GeV,
                  *VertEvent->w_GeV);
}

double SigmaCalc::sigma_tt()
{
  if (*VertEvent->t_GeV < -1.5)
    return 1e-7;
  else{
    double sig_tt = MySigmaTT(*VertEvent->qsq_GeV,
                              -*VertEvent->t_GeV,
                              *VertEvent->w_GeV);
    if ((sig_tt<0) && (Abs(sig_tt)>abs(this->sigma_t())))
      sig_tt = correctionToSigTT(sig_tt,
                                 this->sigma_t(),
                                 *VertEvent->qsq_GeV);
    return sig_tt;
  }
  double sig_tt = MySigmaTT(*VertEvent->qsq_GeV,
                            -*VertEvent->t_GeV,
                            *VertEvent->w_GeV);
  if ((sig_tt<0) && (Abs(sig_tt)>abs(this->sigma_t())))
    sig_tt = correctionToSigTT(sig_tt,
                               this->sigma_t(),
                               *VertEvent->qsq_GeV);
  return sig_tt;
}

double SigmaCalc::sigma_lt()
{
  if (*VertEvent->t_GeV < -1.5)
    return 1e-7;
  else{
    double sig_lt = MySigmaLT(*VertEvent->qsq_GeV,
                              -*VertEvent->t_GeV,
                              *VertEvent->w_GeV);
    if ((sig_lt<0) && (Abs(sig_lt)>abs(this->sigma_t())))
      sig_lt = correctionToSigLT(sig_lt,
                                 this->sigma_t(),
                                 *VertEvent->qsq_GeV);
    return sig_lt;
  }
  double sig_lt = MySigmaLT(*VertEvent->qsq_GeV,
                            -*VertEvent->t_GeV,
                            *VertEvent->w_GeV);
  if ((sig_lt<0) && (Abs(sig_lt)>abs(this->sigma_t())))
    sig_lt = correctionToSigLT(sig_lt,
                               this->sigma_t(),
                               *VertEvent->qsq_GeV);
  return sig_lt;
}

double SigmaCalc::epsilon()
{
  double q = RestEvent->VirtPhot->P()/1000;
  double theta = RestEvent->ScatElec->Theta();

  return 1.0/(1.0 + 2.0*(q*q)/(*VertEvent->qsq_GeV)*
              Power(Tan(theta/2.0),2));

}

double SigmaCalc::sigma_uu()
{

  double phi = *TConEvent->Phi;
  double eps = this->epsilon();
  double siguu = this->sigma_t();
  //cout << siguu << endl;
  siguu += eps*this->sigma_l();
  siguu += Sqrt(2*eps*(1+eps))*this->sigma_lt()*Cos(phi);
  siguu += eps*this->sigma_tt()*Cos(2*phi);
  //siguu /= 2*Pi();
  return siguu;

  //double siguu = ( this->sigma_t() + eps * this->sigma_l() +
  //               eps * TMath::Cos( 2.0 * phi ) * this->sigma_tt() +
  //              TMath::Sqrt( 0.5 * eps * ( 1.0 + eps ) ) * 2.0 * TMath::Cos( phi ) * this->sigma_lt() );
  //return siguu;
}

double SigmaCalc::Sigma_k(int k)
{
  double sigk = this -> Asyms->at(k)
    ->GetAsyAmp(*VertEvent->qsq_GeV,
                *VertEvent->t_prime_GeV);
  // cout << sigk << endl;
  return sigk;
}

double SigmaCalc::sigma_ut()
{
  double phi = *TConEvent->Phi;
  double phi_s = *TConEvent->Phi_s;
  double theta = *RestEvent->Theta;
  double pt = *TConEvent->P_T;

  pt = 0.865; // This is an approximations used in Ahmed's code.
  //It assumes the the q vector is parallel to the lab z axis.

  double sigut = Sin(phi - phi_s)*this->Sigma_k(0);
  //cout << "sigma_ut1\t"<<sigut << endl;
  sigut += Sin(phi_s)*this->Sigma_k(1);
  sigut += Sin(2*phi-phi_s)*this->Sigma_k(2);
  sigut += Sin(phi+phi_s)*this->Sigma_k(3);
  sigut += Sin(3*phi-phi_s)*this->Sigma_k(4);

  //cout << "sigma_ut2\t"<< sigut << endl;

  sigut *= pt/(Sqrt(1-Power(Sin(theta)*Sin(phi_s),2)));
  //cout << "sigma_ut3\t"<< sigut << endl;
  sigut *= this->sigma_uu();

  //cout <<  "sigma_u4\t"<<sigut << endl;


  return sigut;
}

double SigmaCalc::sigma()
{
  double sigma = (this->sigma_uu()+sigma_ut());
  sigma *= this->fluxfactor_col();
  sigma *= this->jacobian_cm() * CofMEvent->ProdPion->P() / (1000*Pi());
  //sigma *= this->jacobian_A();
  sigma *= this->jacobian_cm_col();
  return sigma;
}

double SigmaCalc::weight(int nGen)
{
  return (this->sigma())*PSF()*uBcm2*Lumi/nGen;
}

double SigmaCalc::fluxfactor_col()
{
  double ff = alpha/(2*Power(Pi(),2));
  ff *= (VertEvent->ScatElec->E())/(VertEvent->BeamElec->E());
  ff *= (Power((*VertEvent->w_GeV),2))-Power(neutron_mass_gev,2);
  ff /= 2*neutron_mass_gev;
  ff /= *(VertEvent->qsq_GeV);
  ff /= (1-this->epsilon());
  return ff;
}


double SigmaCalc::jacobian_cm()
{
  double beta = VertEvent->VirtPhot->P()/(VertEvent->VirtPhot->E()+neutron_mass_mev);
  double gamma = (VertEvent->VirtPhot->E() + neutron_mass_mev)/(*VertEvent->w_GeV * 1000);
  double jcm = RestEvent->VirtPhot->P()/1000;
  jcm -= beta * RestEvent->VirtPhot->E()/1000;
  jcm /= gamma * (1-Power(beta,2));
  return jcm;
}

double SigmaCalc::jacobian_cm_col()
{
  double jcol = Power(VertEvent->ProdPion->P(),2)*(*VertEvent->w_GeV*1000);
  jcol /= (CofMEvent->ProdPion->P()
           * Abs((neutron_mass_mev+VertEvent->VirtPhot->E())
                 * VertEvent->ProdPion->P()
                 - (VertEvent->ProdPion->E()*VertEvent->VirtPhot->P()
                    * Cos(VertEvent->ProdPion->Theta())
                    )
                 )
           );
  return jcol;
}

double SigmaCalc::PSF()
{
  double psf = (obj["scat_elec_Emax"].asDouble()-obj["scat_elec_Emin"].asDouble())/1000;
  psf *= (Sin(obj["scat_elec_thetamax"].asDouble()/DEG)
          -Sin(obj["scat_elec_thetamin"].asDouble()/DEG));
  psf *= (Sin(obj["prod_pion_thetamax"].asDouble()/DEG)
          -Sin(obj["prod_pion_thetamin"].asDouble()/DEG));
  psf *= 4*Pi()*Pi();

  return psf;
}

double SigmaCalc::jacobian_A(){
  //return CofMEvent->VirtPhot->Vect().Dot(CofMEvent->ProdPion->Vect())/(1000000*Pi());
  return (this->jacobian_cm() * CofMEvent->ProdPion->P() / (1000*Pi()));
}
