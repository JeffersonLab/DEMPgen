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

class Reaction{

 public:
  Reaction();
  Reaction(TString);
  Reaction(TString, TString);
  ~Reaction();

  void process_reaction();		
  TString GetParticle() {return rParticle;};		
  TString GetHadron() {return rHadron;};
  
 protected:
  TStopwatch tTime;

  TString rParticle;
  TString rHadron;

}; 

	
class PiPlus_Production {

 public:
  PiPlus_Production();
  PiPlus_Production(TString);
  ~PiPlus_Production();

  void process_reaction();		
  TString GetParticle() {return rParticle;};		

 protected:

  void Init();
  void Processing_Event();
  void Progress_Report();
  void Detail_Output();
  void Lund_Output();
  void PiPlus_Pythia6_Out_Init();
  void PiPlus_Pythia6_Output();
  void PiPlus_HEPMC3_Out_Init();
  void PiPlus_HEPMC3_Output();

  TRandom2* rRand;
		
  Particle_t recoil_nucleon;
  Particle_t produced_X;

  Double_t Get_Phi_X_LeptonPlane_RF();
  Double_t Get_Phi_TargPol_LeptonPlane_RF();

  Double_t Get_Total_Cross_Section(); 

  Double_t GetPi0_CrossSection();
  Double_t GetPiPlus_CrossSection();

  /*--------------------------------------------------*/
  // Parameters

  TStopwatch tTime;
  TString rParticle;
  TString rParticle_charge;
  TString rParticle_scat_nucleon;

  std::string sTFile;   /// Generator output files. For documentation and monitoring purposes 
  std::string sLFile;   /// Lund input file into the EIC simulation

  std::ofstream ppiOut;     
  std::ofstream ppiDetails;
		
  long long int qsq_ev, t_ev, w_neg_ev, w_ev;
		
  long long int rNEvents;
  long long int rNEvent_itt;
  TDatime dFractTime; 

  double rDEG2RAD;
                   
  double fX_Theta_I, fX_Theta_F;

  TLorentzVector GetProtonVector_lab();
  TLorentzVector GetElectronVector_lab();

  TLorentzVector r_lproton;     // Proton in collider (lab) frame
  TLorentzVector r_lprotong;	

  TLorentzVector r_lelectron;   // Electron in collider (lab) frame
  TLorentzVector r_lelectrong;	
		
  TVector3 beta_col_rf; // Boost vector from collider (lab) frame to protons rest frame (Fix target)

  void Consider_Proton_Fermi_Momentum();

  Double_t rFermiMomentum;

  Double_t fX_Theta_Col, fX_Phi_Col;

  TLorentzVector r_lscatelec; 
  TLorentzVector r_lscatelecg;

  TLorentzVector r_lphoton;
  TLorentzVector r_lphotong;
    
  TLorentzVector r_lX;
  TLorentzVector r_lX_g;

  double fX_Mass;
  double fX_Mass_GeV;

  double f_Scat_Nucleon_Mass;     
  double f_Scat_Nucleon_Mass_GeV;

  TLorentzVector r_l_scat_nucleon;
  TLorentzVector r_l_scat_nucleon_g;

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

  TLorentzVector lX_rf;
  TLorentzVector lX_rfg;

  TLorentzVector l_scat_nucleon_rf;
  TLorentzVector l_scat_nucleon_rf_g;

  ///////////////////////////////////////////
  /// Center of Mass parameters for particle X

  double fBeta_CM_RF, fGamma_CM_RF, fX_Energy_CM, fX_Mom_CM, fX_Energy_CM_GeV, fX_Mom_CM_GeV;

  TLorentzVector lt;
  TLorentzVector ltg;
		
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

}; 

//	
//class PiPlus_Production: public virtual Reaction{
//
//	public:
//		PiPlus_Production();
//		PiPlus_Production(TString);
//		~PiPlus_Production();
//
//		void process_reaction();		
// 		void Lund_Output();
//
//     protected:
//
// 		void Init();
// 		void Processing_Event();
// 		void Progress_Report();
//		void Detail_Output();
//
//
//
// 
//// }
//

class Pi0_Production:PiPlus_Production{
 
 public:
  Pi0_Production();
  Pi0_Production(TString);
  ~Pi0_Production();

  void process_reaction();
  void Detail_Output();
  void Processing_Event();
  Double_t Get_CrossSection();

  void Pi0_decay(TLorentzVector);

  bool if_pi0_decay;

  void Pi0_Decay_Pythia6_Out_Init();
  void Pi0_Decay_Pythia6_Output();

  unsigned long long int print_itt;

 private:

  Double_t theta_X_rf;

  Double_t ft_min;
  Double_t fu_min;

  Double_t ft;
  Double_t fu;

  std::ofstream polar_out;     

  TLorentzVector l_photon_1;
  TLorentzVector l_photon_2;

  void Pi0_Lund_Output();
  void Pi0_Decay_Lund_Output();

  //  		template <class T> 
  //  		inline int 
  //  		sgn(T val) {
  //      		return (T(0) < val) - (val < T(0));
  //  		}

  template <class T>
    T Sign (T a, T b) {
    return (a < b) - (b < a);
  }

};

class KPlus_Production {
  
 public:
                 
  KPlus_Production();
  KPlus_Production(TString, TString);
  ~KPlus_Production();

  void process_reaction();
  TString GetParticle() {return rParticle;};
  TString GetHadron() {return rHadron;};

 protected:

  void Init();
  void Processing_Event();
  void Progress_Report();
  void Detail_Output();
  void Lund_Output();
  void KPlus_Pythia6_Out_Init();
  void KPlus_Pythia6_Output();
  void KPlus_HEPMC3_Out_Init();
  void KPlus_HEPMC3_Output();
  TRandom2* rRand;
		  
  Particle_t recoil_hadron;                  
  Particle_t produced_X;

  Double_t Get_Phi_X_LeptonPlane_RF();
  Double_t Get_Phi_TargPol_LeptonPlane_RF();
  Double_t Get_Total_Cross_Section();
  Double_t GetPi0_CrossSection();
  Double_t GetPiPlus_CrossSection();
  Double_t GetKPlus_CrossSection();
  // SJDK - 08/02/22 - Only one cross section value for now? Split into KLambda/KSigma later
  //Double_t GetKLambda_CrossSection();
  //Double_t GetKSigma_CrossSection();

  /*--------------------------------------------------*/
  // Parameters

  TStopwatch tTime;
  TString rParticle;
  TString rParticle_charge;
  TString rParticle_scat_hadron;
  TString rHadron;

  std::string sTFile;   /// Generator output files. For documentation and monitoring purposes
  std::string sLFile;   /// Lund input file into the EIC simulation
  std::ofstream ppiOut;
  std::ofstream ppiDetails;

  int qsq_ev, t_ev, w_neg_ev, w_ev;
  long long int rNEvents;
  long long int rNEvent_itt;
  TDatime dFractTime;
  double rDEG2RAD;

  double fX_Theta_I, fX_Theta_F;

  TLorentzVector GetProtonVector_lab();
  TLorentzVector GetElectronVector_lab();
  TLorentzVector r_lproton;     // Proton in collider (lab) frame
  TLorentzVector r_lprotong;
  TLorentzVector r_lelectron;   // Electron in collider (lab) frame
  TLorentzVector r_lelectrong;
	  
  TVector3 beta_col_rf; // Boost vector from collider (lab) frame to protons rest frame (Fix target)
  void Consider_Proton_Fermi_Momentum();
  Double_t rFermiMomentum;
  Double_t fX_Theta_Col, fX_Phi_Col;
  TLorentzVector r_lscatelec;
  TLorentzVector r_lscatelecg;
  TLorentzVector r_lphoton;
  TLorentzVector r_lphotong;
  TLorentzVector r_lX;
  TLorentzVector r_lX_g;
  double fX_Mass;
  double fX_Mass_GeV;
  double f_Scat_hadron_Mass;
  double f_Scat_hadron_Mass_GeV;

  TLorentzVector r_l_scat_hadron;
  TLorentzVector r_l_scat_hadron_g;
  TLorentzVector r_lw;
  TLorentzVector lwp;
  TLorentzVector fsini;
  TLorentzVector fsfin;
  TLorentzVector fsinig;
  TLorentzVector fsfing;
  pim* pd;
	  
  ///////////////////////////////////////////
  // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss
                       
  TLorentzVector lproton_rf;
  TLorentzVector lproton_rfg;
  TLorentzVector lelectron_rf;
  TLorentzVector lelectron_rfg;
  TLorentzVector lscatelec_rf;
  TLorentzVector lscatelec_rfg;
  TLorentzVector lphoton_rf;
  TLorentzVector lphoton_rfg;
  TLorentzVector lX_rf;
  TLorentzVector lX_rfg;
  TLorentzVector l_scat_hadron_rf;
  TLorentzVector l_scat_hadron_rf_g;

  ///////////////////////////////////////////
  /// Center of Mass parameters for particle X

  double fBeta_CM_RF, fGamma_CM_RF, fX_Energy_CM, fX_Mom_CM, fX_Energy_CM_GeV, fX_Mom_CM_GeV;
	  
  TLorentzVector lt;
  TLorentzVector ltg;
      
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

};

class DEMP_Reaction {

 public:
  DEMP_Reaction();
  DEMP_Reaction(TString, TString);
  ~DEMP_Reaction();

  void process_reaction();		
  TString GetParticle() {return rParticle;};		
  TString GetHadron() {return rHadron;};

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
		
  Particle_t recoil_hadron;
  Particle_t produced_X;

  Double_t Get_Phi_X_LeptonPlane_RF();
  Double_t Get_Phi_TargPol_LeptonPlane_RF();

  Double_t Get_Total_Cross_Section(); 

  Double_t GetPi0_CrossSection();
  Double_t GetPiPlus_CrossSection();
  Double_t GetKPlus_CrossSection();


  /*--------------------------------------------------*/
  // Parameters

  TStopwatch tTime;
  TString rParticle;
  TString rParticle_charge;
  TString rParticle_scat_hadron;
  TString rHadron;

  std::string sTFile;   /// Generator output files. For documentation and monitoring purposes 
  std::string sLFile;   /// Lund input file into the EIC simulation

  std::ofstream DEMPOut;     
  std::ofstream DEMPDetails;
		
  long long int qsq_ev, t_ev, w_neg_ev, w_ev;
		
  long long int rNEvents;
  long long int rNEvent_itt;
  TDatime dFractTime; 

  double rDEG2RAD;
                   
  double fX_Theta_I, fX_Theta_F;

  TLorentzVector GetProtonVector_lab();
  TLorentzVector GetElectronVector_lab();

  TLorentzVector r_lproton;     // Proton in collider (lab) frame
  TLorentzVector r_lprotong;	

  TLorentzVector r_lelectron;   // Electron in collider (lab) frame
  TLorentzVector r_lelectrong;	
		
  TVector3 beta_col_rf; // Boost vector from collider (lab) frame to protons rest frame (Fix target)

  void Consider_Proton_Fermi_Momentum();

  Double_t rFermiMomentum;

  Double_t fX_Theta_Col, fX_Phi_Col;

  TLorentzVector r_lscatelec; 
  TLorentzVector r_lscatelecg;

  TLorentzVector r_lphoton;
  TLorentzVector r_lphotong;
    
  TLorentzVector r_lX;
  TLorentzVector r_lX_g;

  double fX_Mass;
  double fX_Mass_GeV;

  double f_Scat_hadron_Mass;     
  double f_Scat_hadron_Mass_GeV;

  TLorentzVector r_l_scat_hadron;
  TLorentzVector r_l_scat_hadron_g;

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

  TLorentzVector lX_rf;
  TLorentzVector lX_rfg;

  TLorentzVector l_scat_hadron_rf;
  TLorentzVector l_scat_hadron_rf_g;

  ///////////////////////////////////////////
  /// Center of Mass parameters for particle X

  double fBeta_CM_RF, fGamma_CM_RF, fX_Energy_CM, fX_Mom_CM, fX_Energy_CM_GeV, fX_Mom_CM_GeV;

  TLorentzVector lt;
  TLorentzVector ltg;
		
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

};

# endif
