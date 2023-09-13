# ifndef REACTION_ROUNTINE_CC
# define REACTION_ROUNTINE_CC

#include "eic_pim.h"

#include <string>
#include <fstream>

#include <TStopwatch.h>
#include <TDatime.h>

#include "TF1.h"
#include "TLorentzVector.h"

#include "particleType.h"

#include "TCanvas.h"

#include "Particle.hxx"
#include "CustomRand.hxx"


class Reaction{

 public:
  Reaction();
  Reaction(TString);
  Reaction(TString, TString);
  ~Reaction();

  void process_reaction();		
  TString GetEjectile() {return rEjectile;};		
  TString GetRecoilHadron() {return rRecoil;};
  
 protected:
  TStopwatch tTime;

  TString rEjectile;
  TString rRecoil;

};

class DEMP_Reaction {

 public:
  DEMP_Reaction();
  DEMP_Reaction(TString, TString);
  ~DEMP_Reaction();

  void process_reaction();		
  TString GetEjectile() {return rEjectile;};		
  TString GetRecoilHadron() {return rRecoil;};

 protected:

  void Init();
  void Processing_Event();
  void Progress_Report();
  void Detail_Output();
  void Lund_Output();
  void DEMPReact_Pythia6_Out_Init();
  void DEMPReact_Pythia6_Output();
  void DEMPReact_HEPMC3_Out_Init();
  void DEMPReact_HEPMC3_Output();

  TRandom2* rRanBd;
  TRandom2* rRand;
		
  Particle_t recoil_hadron;
  Particle_t produced_X;

  Double_t Get_Phi_X_LeptonPlane_RF();
  Double_t Get_Phi_TargPol_LeptonPlane_RF();

  Double_t Get_Total_Cross_Section(); 

  //  Double_t GetPi0_CrossSection();

  /*--------------------------------------------------*/
  // Parameters

  TStopwatch tTime;
  TString rEjectile;
  TString rEjectile_charge;
  TString rEjectile_scat_hadron;
  TString rRecoil;

  std::string sTFile;   /// Generator output files. For documentation and monitoring purposes 
  std::string sLFile;   /// Lund input file into the EIC simulation
  std::string sDFile;   /// Root diagnostic plot in root file format

  std::ofstream DEMPOut;     
  std::ofstream DEMPDetails;

  TFile* dRootFile;
  TTree* dRootTree;
		
  long long int qsq_ev, t_ev, w_neg_ev, w_ev;
		
  long long int rNEvents;
  long long int rNEvent_itt;
  TDatime dFractTime; 

  double rDEG2RAD;
                   
  double f_Ejectile_Theta_I, f_Ejectile_Theta_F;

  TLorentzVector GetProtonVector_lab();
  TLorentzVector GetElectronVector_lab();

  TLorentzVector r_lproton;     // Proton in collider (lab) frame
  TLorentzVector r_lprotong;	

  TLorentzVector r_lelectron;   // Electron in collider (lab) frame
  TLorentzVector r_lelectrong;	

  TLorentzVector r_lhadron_beam;	
  TLorentzVector r_lhadron_beamg;

  Double_t r_lhadron_beam_mass;	
  Double_t r_lhadron_beam_massg;	
		
  TVector3 beta_col_rf; // Boost vector from collider (lab) frame to protons rest frame (Fix target)

  void Consider_Proton_Fermi_Momentum();

  Double_t rFermiMomentum;

  Double_t f_Ejectile_Theta_Col, f_Ejectile_Phi_Col;

  TLorentzVector r_lscatelec; 
  TLorentzVector r_lscatelecg;

  TLorentzVector r_lphoton;
  TLorentzVector r_lphotong;
    
  TLorentzVector r_l_Ejectile;
  TLorentzVector r_l_Ejectile_g;

  TLorentzVector r_l_Ejectile_solved;
  TLorentzVector r_l_Recoil_solved;

  double f_Ejectile_Mass;
  double f_Ejectile_Mass_GeV;

  double f_Recoil_Mass;     
  double f_Recoil_Mass_GeV;

  TLorentzVector l_Recoil;
  TLorentzVector l_Recoil_g;

  TLorentzVector r_lw;

  TLorentzVector lwp;

  TLorentzVector fsini;
  TLorentzVector fsfin;
  TLorentzVector fsinig;
  TLorentzVector fsfing;

  pim* pd;

  ///////////////////////////////////////////
  //Transformation of e', pi- and recoil proton to target's rest frmae without energy loss 

  TLorentzVector lproton_rf;
  TLorentzVector lproton_rfg;

  TLorentzVector lelectron_rf;
  TLorentzVector lelectron_rfg;

  TLorentzVector lscatelec_rf;
  TLorentzVector lscatelec_rfg;

  TLorentzVector lphoton_rf;
  TLorentzVector lphoton_rfg;

  TLorentzVector l_Ejectile_rf;
  TLorentzVector l_Ejectile_rfg;

  TLorentzVector l_scat_hadron_rf;
  TLorentzVector l_scat_hadron_rf_g;

  ///////////////////////////////////////////
  /// Center of Mass parameters for particle X

  double fBeta_CM_RF, fGamma_CM_RF, f_Ejectile_Energy_CM, f_Ejectile_Mom_CM, f_Ejectile_Energy_CM_GeV, f_Ejectile_Mom_CM_GeV;

  TLorentzVector lt;
  TLorentzVector ltg;
  TLorentzVector lu;
  TLorentzVector lug;

  ///////////////////////////////////////////

  TVector3 v3Photon;   
  TVector3 v3Electron; 
  TVector3 v3X;    
  TVector3 v3S;        
  TVector3 v3PhotonUnit;
  TVector3 v3QxL;       
  TVector3 v3QxP;       
  TVector3 v3QxS;       
  TVector3 v3LxP;       
  TVector3 v3LxS;       
  TVector3 v3PxL;       
  TVector3 v3QUnitxL;  
  TVector3 v3QUnitxP;   
  TVector3 v3QUnitxS;   

  double fCos_Phi_X_LeptonPlane_RF, fSin_Phi_X_LeptonPlane_RF, fTheta_X_Photon_RF, fPhi_X_LeptonPlane_RF;
  Double_t r_fSig;
  Double_t r_fSig_T;
  Double_t r_fSig_L;

  unsigned long long int print_itt;

  ///*--------------------------------------------------*/ 
  // Rory Check algorithm

  TLorentzVector* Interaction_Solve;
  TLorentzVector* Target_Solve;

  TLorentzVector* VertBeamElec;
  TLorentzVector* VertScatElec;
  
  TLorentzVector* Initial;
  TLorentzVector* Target;
  TLorentzVector* Photon;
  TLorentzVector* Interaction;
  TLorentzVector* Final;
  
  bool SolnCheck();
  double W_in_Solve(); 
  double W_out_Solve();
  double W_in_val;

  TRandom3* CoinToss;
  CustomRand* AngleGen;
 
  TF1* F;
  TVector3* UnitVect;

  int Solve();
  int Solve(double theta, double phi);

  double W_in();
  double W_out();
  
  ///*--------------------------------------------------*/ 
  // Needed for the Solve function

  double theta;
  double phi;
  double P;
  double P2;

  double tc;
  double tc_GeV;
  double uc;
  double uc_GeV;
  ///*--------------------------------------------------*/ 

};

//class Pi0_Production:DEMP_Reaction{
// 
// public:
//  Pi0_Production();
//  Pi0_Production(TString);
//  ~Pi0_Production();
//
//  void process_reaction();
//  void Detail_Output();
//  void Processing_Event();
//  Double_t Get_CrossSection();
//
//  void Pi0_decay(TLorentzVector);
//
//  bool if_pi0_decay;
//
//  void Pi0_Decay_Pythia6_Out_Init();
//  void Pi0_Decay_Pythia6_Output();
//  
//  ///----------------------------------------------------*/
//  /// Output algorithm into HEPMC3 format
//
//  void Pi0_HEPMC3_Out_Init();
//  void Pi0_HEPMC3_Output();
//
//  unsigned long long int print_itt;
//
// private:
//
//  Double_t theta_X_rf;
//
//  Double_t ft_min;
//  Double_t fu_min;
//
//  Double_t ft;
//  Double_t fu;
//
//  std::ofstream polar_out;     
//
//  TLorentzVector l_photon_1;
//  TLorentzVector l_photon_2;
//
//  void Pi0_Lund_Output();
//  void Pi0_Decay_Lund_Output();
//
//  //  		template <class T> 
//  //  		inline int 
//  //  		sgn(T val) {
//  //      		return (T(0) < val) - (val < T(0));
//  //  		}
//
//  template <class T>
//    T Sign (T a, T b) {
//    return (a < b) - (b < a);
//  }
//
//};

# endif
