# ifndef EIC_PIM_H
# define EIC_PIM_H

#include <iostream>
#include <string>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"

using std::vector;

class pim {

public:
  pim(); 
  pim(int);

  void Initilize();
  int CheckLaws(TLorentzVector P_E0, TLorentzVector P_t, TLorentzVector P_e, TLorentzVector P_pim, TLorentzVector P_pro);
  int CheckLaws(TLorentzVector P_E0, TLorentzVector P_t, TLorentzVector P_e, TLorentzVector P_pim, TLorentzVector P_pro, double fDiff_E);
  void setrootfile(std::string myRootFile );
  double fermiMomentum();

private:
  Int_t gen_seed = 0;
	
  std::string pParticle;
  std::string pcharge;

  /* double correctedPhi(); */
  /* double correctedPhiS(); */

};

//extern TRandom2 *fRandom;                    

extern TRandom3 *fRandom;                    

extern TFile *f;

extern TTree *t1;

extern int gKinematics_type;
extern TString gfile_name;
extern TString gParticle;
extern TString gHadron;
extern bool gPi0_decay;
extern bool UseSolve;
extern std::string gDet_location;
extern std::string gOutputType;
extern std::string gBeamPart;
extern float fProton_incidence_phi;

extern int fSeed;

extern bool allset;
extern bool kCalcFermi;
extern bool kCalcBremss;
extern bool kCalcIon;
extern bool kCalcBremssEle;
extern bool kCalcIonEle;
extern bool kSConserve;
extern bool kFSI;
extern bool kMSele;
extern bool kMS;

extern double fKaon_Mass;
extern double fKaon_Mass_GeV;

extern double fLambda_Mass;                             
extern double fLambda_Mass_GeV;

extern double fSigma_Mass;
extern double fSigma_Mass_GeV;

extern double fOmega_Mass; 
extern double fOmega_Mass_GeV; 

extern double fOmega_Theta_Col; 
extern double fOmega_Phi_Col; 

extern double fOmega_Theta_I; 
extern double fOmega_Theta_F; 

extern double fOmega_Energy_CM;    
extern double fOmega_Mom_CM;       
extern double fOmega_Energy_CM_GeV;
extern double fOmega_Mom_CM_GeV;   

extern double fPhi_Omega_LeptonPlane_RF;
extern double fCos_Phi_Omega_LeptonPlane_RF; 
extern double fSin_Phi_Omega_LeptonPlane_RF;
extern double fTheta_Omega_Photon_RF;

extern int fWLessShell;
extern int fWLess1P9;
extern int fSDiff;

//extern long int fNEvents;

extern unsigned long long int fNEvents;
extern unsigned long long int fNRecorded;
extern unsigned long long int fNGenerated;
extern unsigned long long int fWSqNeg;
extern unsigned long long int fNMomConserve;
extern unsigned long long int fNSigmaNeg;
extern unsigned long long int fNWeightNeg;

extern unsigned long long int fNaN;
extern unsigned long long int fConserve;

extern unsigned long long int fNWeightUnphys;
extern unsigned long long int fNWeightReject;

extern unsigned long long int fSolveEvents_0Sol;
extern unsigned long long int fSolveEvents_1Sol;
extern unsigned long long int fSolveEvents_2Sol;

extern unsigned long long int fNFile;

extern double fK;
extern double fm;
extern double fElectron_Kin_Col_GeV;
extern double fElectron_Kin_Col;
extern double fRand;
extern double fLumi;
extern double fuBcm2;
extern double fPI;
extern double fDEG2RAD;
extern double fRAD2DEG;
extern double fEBeam;
extern double fPBeam;
// 13/09/23 - SJDK - New generic HBeam value (rather than proton beam)
extern double fHBeam;
extern double fScatElec_Theta_I;
extern double fScatElec_Theta_F;
extern double fPion_Theta_I; // SJDK 19/12/22 - These should be removed in future, specific to pion reaction cases. Should be generic MesonX
extern double fPion_Theta_F;
extern double fEjectileX_Theta_I;
extern double fEjectileX_Theta_F;
extern double fScatElec_E_Hi;
extern double fScatElec_E_Lo;
extern double fPSF;
// SJDK 03/04/23 - New Variables for min/max allowed Q2/W, set by particle type - Add an override in .json file?
extern double fQsq_Min;
extern double fQsq_Max;
extern double fW_Min;
extern double fW_Max;
extern double fT_Max;

extern double fMandSConserve;
extern double fTop_Pion_Mom;
extern double fBot_Pion_Mom;
extern double fPion_Mom_Same;
extern double fEnergyConserve;
extern double fXMomConserve;
extern double fYMomConserve;
extern double fZMomConserve;
extern double fXMomConserve_RF;
extern double fYMomConserve_RF;
extern double fZMomConserve_RF;
extern double fEnergyConserve_RF;

extern double fDiff;
extern double fRatio;
extern double fPion_Alpha;
extern double fPion_Beta;
extern double fS_I_RF;
extern double fS_F_RF;
extern double fS_I_Col;
extern double fS_F_Col;
extern double fS_I_RF_GeV;
extern double fS_F_RF_GeV;
extern double fS_I_Col_GeV;
extern double fS_F_Col_GeV;

extern double fProton_Energy_Col;
extern double fProton_Mom_Col;
extern double fProton_Theta_Col;
extern double fProton_Phi_Col;
extern double fProton_MomZ_Col;
extern double fProton_MomX_Col;
extern double fProton_MomY_Col;
extern double fProton_Energy_Col_GeV;
extern double fProton_Mom_Col_GeV;
extern double fProton_MomX_Col_GeV;
extern double fProton_MomY_Col_GeV;
extern double fProton_MomZ_Col_GeV;

extern double fFSIProton_Energy_Col;
extern double fFSIProton_Mom_Col;
extern double fFSIProton_Theta_Col;
extern double fFSIProton_Phi_Col;
extern double fFSIProton_MomZ_Col;
extern double fFSIProton_MomX_Col;
extern double fFSIProton_MomY_Col;
extern double fFSIProton_Energy_Col_GeV;
extern double fFSIProton_Mom_Col_GeV;
extern double fFSIProton_MomX_Col_GeV;
extern double fFSIProton_MomY_Col_GeV;
extern double fFSIProton_MomZ_Col_GeV;

extern double fTarget_Energy_Col;
extern double fTarget_Mom_Col;
extern double fTarget_Theta_Col;
extern double fTarget_Phi_Col;
extern double fTarget_MomZ_Col;
extern double fTarget_MomX_Col;
extern double fTarget_MomY_Col;
extern double fTarget_Energy_Col_GeV;
extern double fTarget_Mom_Col_GeV;
extern double fTarget_MomX_Col_GeV;
extern double fTarget_MomY_Col_GeV;
extern double fTarget_MomZ_Col_GeV;

extern double fTarget_Pol0_Col;
extern double fTarget_PolX_Col;
extern double fTarget_PolY_Col;
extern double fTarget_PolZ_Col;
extern double fTarget_Pol0_RF;
extern double fTarget_PolX_RF;
extern double fTarget_PolY_RF;
extern double fTarget_PolZ_RF;

extern double fBetaX_Col_RF;
extern double fBetaY_Col_RF;
extern double fBetaZ_Col_RF;
extern double fBeta_Col_RF;
extern double fGamma_Col_RF;

extern double fProton_MomX_RF;
extern double fProton_MomY_RF;
extern double fProton_MomZ_RF;
extern double fProton_Mom_RF;
extern double fProton_Energy_RF;
extern double fProton_MomX_RF_GeV;
extern double fProton_MomY_RF_GeV;
extern double fProton_MomZ_RF_GeV;
extern double fProton_Mom_RF_GeV;
extern double fProton_Energy_RF_GeV;

extern double fScatElec_Angle;    
extern double fScatElec_Alpha_RF;    
extern double fScatElec_Beta_RF;

extern double fVertex_X;
extern double fVertex_Y;
extern double fVertex_Z;
extern double fProton_Kin_Col_GeV;
extern double fElectron_Mass;
extern double fElectron_Mass_GeV;
extern double fProton_Mass;
extern double fProton_Mass_GeV;
extern double fNeutron_Mass;
extern double fNeutron_Mass_GeV;
extern double fPion_Mass;
extern double fPion_Mass_GeV;
extern double fPiion_Phi;
extern double fAlpha;
extern double fPi;
extern double fMom_Ratio;
extern double fMom_Dif;
extern double fPionEnergyCMLess;
extern double fSNotEqual;
extern double fMode_Epsi;
extern double fRecoilProton_Mass;
extern double fRecoilProton_Mass_GeV;

extern double fElectron_Energy_Col;
extern double fElectron_MomZ_Col;
extern double fElectron_MomX_Col;
extern double fElectron_MomY_Col;
extern double fElectron_Theta_Col;
extern double fElectron_Phi_Col;
extern double fElectron_Mom_Col;

extern double fElectron_MS_Energy_Col;
extern double fElectron_MS_MomZ_Col;
extern double fElectron_MS_MomX_Col;
extern double fElectron_MS_MomY_Col;
extern double fElectron_MS_Theta_Col;
extern double fElectron_MS_Phi_Col;
extern double fElectron_MS_Mom_Col;

extern double fElectron_Energy_Col_GeV;
extern double fElectron_Mom_Col_GeV;
extern double fElectron_MomX_Col_GeV;
extern double fElectron_MomY_Col_GeV;
extern double fElectron_MomZ_Col_GeV;
extern double fElectronEnergyLess;
extern double fElectronThetaLess;
extern double fRadiation_Lenght_Air;

extern double fElectron_Targ_Thickness;
extern double fElectron_Targ_Thickness_RadLen;
extern double fElectron_Targ_BT;
extern double fElectron_Targ_Bremss_Loss;
extern double fElectron_Targ_Ion_Loss;
extern double fElectron_TargWindow_Bremss_Loss;
extern double fElectron_TargWindow_Ion_Loss;

extern double fElectron_Air_Thickness;
extern double fElectron_Air_Thickness_RadLen;
extern double fElectron_Air_BT;
extern double fElectron_Air_Bremss_Loss;
extern double fElectron_Air_Ion_Loss;
extern double fElectron_Corrected_Energy_Col;
extern double fElectron_Corrected_Mom_Col;
extern double fElectron_Corrected_MomX_Col;
extern double fElectron_Corrected_MomY_Col;
extern double fElectron_Corrected_MomZ_Col;
extern double fElectron_Corrected_Theta_Col;
extern double fElectron_Corrected_Phi_Col;
extern double fElectron_Delta_Mom_Col;
extern double fElectron_Corrected_Energy_Col_GeV;
extern double fElectron_Corrected_Mom_Col_GeV;
extern double fElectron_Corrected_MomX_Col_GeV;
extern double fElectron_Corrected_MomY_Col_GeV;
extern double fElectron_Corrected_MomZ_Col_GeV;
extern double fElectron_Delta_Mom_Col_GeV;

extern double fScatElec_MS_Energy_Col;
extern double fScatElec_MS_MomZ_Col;
extern double fScatElec_MS_MomX_Col;
extern double fScatElec_MS_MomY_Col;
extern double fScatElec_MS_Theta_Col;
extern double fScatElec_MS_Phi_Col;
extern double fScatElec_MS_Mom_Col;

extern double fScatElec_Energy_Col;
extern double fScatElec_MomZ_Col;
extern double fScatElec_MomX_Col;
extern double fScatElec_MomY_Col;
extern double fScatElec_Theta_Col;
extern double fScatElec_Phi_Col;
extern double fScatElec_Mom_Col;
extern double fScatElec_Energy_Col_GeV;
extern double fScatElec_Mom_Col_GeV;
extern double fScatElec_MomX_Col_GeV;
extern double fScatElec_MomY_Col_GeV;
extern double fScatElec_MomZ_Col_GeV;
extern double fScatElecEnergyLess;
extern double fScatElecThetaLess;
extern double fScatElec_Targ_Thickness;
extern double fScatElec_Targ_Thickness_RadLen;
extern double fScatElec_Targ_BT;
extern double fScatElec_Targ_Bremss_Loss;
extern double fScatElec_Targ_Ion_Loss;
extern double fScatElec_Air_Thickness;
extern double fScatElec_Air_Thickness_RadLen;
extern double fScatElec_Air_BT;
extern double fScatElec_Air_Bremss_Loss;
extern double fScatElec_Air_Ion_Loss;
extern double fScatElec_Corrected_Energy_Col;
extern double fScatElec_Corrected_Mom_Col;
extern double fScatElec_Corrected_MomX_Col;
extern double fScatElec_Corrected_MomY_Col;
extern double fScatElec_Corrected_MomZ_Col;
extern double fScatElec_Corrected_Theta_Col;
extern double fScatElec_Corrected_Phi_Col;
extern double fScatElec_Delta_Mom_Col;
extern double fScatElec_Corrected_Energy_Col_GeV;
extern double fScatElec_Corrected_Mom_Col_GeV;
extern double fScatElec_Corrected_MomX_Col_GeV;
extern double fScatElec_Corrected_MomY_Col_GeV;
extern double fScatElec_Corrected_MomZ_Col_GeV;
extern double fScatElec_Delta_Mom_Col_GeV;
extern double fScatElec_TargWindow_Bremss_Loss;
extern double fScatElec_TargWindow_Ion_Loss;
extern double fTargWindow_Thickness;
extern double fTargWindow_Thickness_RadLen;
extern double fTargWindow_BT;

extern double fPion_TargWindow_Ion_Loss;
extern double fPion_Targ_Thickness;
extern double fPion_Targ_Thickness_RadLen;
extern double fPion_Targ_BT;
extern double fPion_Targ_Bremss_Loss;
extern double fPion_Targ_Ion_Loss;
extern double fPion_Air_Thickness;
extern double fPion_Air_Thickness_RadLen;
extern double fPion_Air_BT;
extern double fPion_Air_Bremss_Loss;
extern double fPion_Air_Ion_Loss;

extern double fPion_MS_Energy_Col;
extern double fPion_MS_MomZ_Col;
extern double fPion_MS_MomX_Col;
extern double fPion_MS_MomY_Col;
extern double fPion_MS_Theta_Col;
extern double fPion_MS_Phi_Col;
extern double fPion_MS_Mom_Col;

extern double fPion_Theta_Col;
extern double fPion_Phi_Col;
extern double fPion_Energy_Col;
extern double fPion_Mom_Col;
extern double fPion_MomZ_Col;
extern double fPion_MomX_Col;
extern double fPion_MomY_Col;
extern double fPion_Energy_Col_GeV;
extern double fPion_Mom_Col_GeV;
extern double fPion_MomX_Col_GeV;
extern double fPion_MomY_Col_GeV;
extern double fPion_MomZ_Col_GeV;

extern double fPion_FSI_Theta_Col;
extern double fPion_FSI_Phi_Col;
extern double fPion_FSI_Energy_Col;
extern double fPion_FSI_Mom_Col;
extern double fPion_FSI_MomZ_Col;
extern double fPion_FSI_MomX_Col;
extern double fPion_FSI_MomY_Col;
extern double fPion_FSI_Energy_Col_GeV;
extern double fPion_FSI_Mom_Col_GeV;
extern double fPion_FSI_MomX_Col_GeV;
extern double fPion_FSI_MomY_Col_GeV;
extern double fPion_FSI_MomZ_Col_GeV;

extern double fPion_Corrected_Theta_Col;
extern double fPion_Corrected_Phi_Col;
extern double fPion_Corrected_Energy_Col;
extern double fPion_Corrected_Mom_Col;
extern double fPion_Corrected_MomX_Col;
extern double fPion_Corrected_MomY_Col;
extern double fPion_Corrected_MomZ_Col;
extern double fPion_Delta_Mom_Col;
extern double fPion_Corrected_Energy_Col_GeV;
extern double fPion_Corrected_Mom_Col_GeV;
extern double fPion_Corrected_MomX_Col_GeV;
extern double fPion_Corrected_MomY_Col_GeV;
extern double fPion_Corrected_MomZ_Col_GeV;
extern double fPion_Delta_Mom_Col_GeV;

extern double fKaon_Theta_Col; 
extern double fKaon_Phi_Col;
extern double fKaon_Energy_Col;
extern double fKaon_Mom_Col;
extern double fKaon_MomZ_Col;
extern double fKaon_MomX_Col;    
extern double fKaon_MomY_Col;
extern double fKaon_Energy_Col_GeV;
extern double fKaon_Mom_Col_GeV;
extern double fKaon_MomX_Col_GeV;
extern double fKaon_MomY_Col_GeV;   
extern double fKaon_MomZ_Col_GeV;
extern double fScathad_Theta_Col;
extern double fScathad_Phi_Col;
extern double fScathad_Energy_Col;
extern double fScathad_Mom_Col;
extern double fScathad_MomZ_Col;
extern double fScathad_MomX_Col;
extern double fScathad_MomY_Col;
extern double fScathad_Energy_Col_GeV;
extern double fScathad_Mom_Col_GeV;
extern double fScathad_MomX_Col_GeV;  
extern double fScathad_MomY_Col_GeV; 
extern double fScathad_MomZ_Col_GeV;

extern double fNeutron_MS_Energy_Col;
extern double fNeutron_MS_MomZ_Col;
extern double fNeutron_MS_MomX_Col;
extern double fNeutron_MS_MomY_Col;
extern double fNeutron_MS_Theta_Col;
extern double fNeutron_MS_Phi_Col;
extern double fNeutron_MS_Mom_Col;

extern double fNeutron_TargWindow_Ion_Loss;
extern double fNeutron_Targ_Thickness;
extern double fNeutron_Targ_Thickness_RadLen;
extern double fNeutron_Targ_BT;
extern double fNeutron_Targ_Bremss_Loss;
extern double fNeutron_Targ_Ion_Loss;
extern double fNeutron_Air_Thickness;
extern double fNeutron_Air_Thickness_RadLen;
extern double fNeutron_Air_BT;
extern double fNeutron_Air_Bremss_Loss;
extern double fNeutron_Air_Ion_Loss;
extern double fNeutron_Theta_Col;
extern double fNeutron_Phi_Col;
extern double fNeutron_Energy_Col;
extern double fNeutron_Mom_Col;
extern double fNeutron_MomZ_Col;
extern double fNeutron_MomX_Col;
extern double fNeutron_MomY_Col;
extern double fNeutron_Energy_Col_GeV;
extern double fNeutron_Mom_Col_GeV;
extern double fNeutron_MomX_Col_GeV;
extern double fNeutron_MomY_Col_GeV;
extern double fNeutron_MomZ_Col_GeV;
extern double fNeutron_Corrected_Theta_Col;
extern double fNeutron_Corrected_Phi_Col;
extern double fNeutron_Corrected_Energy_Col;
extern double fNeutron_Corrected_Mom_Col;
extern double fNeutron_Corrected_MomX_Col;
extern double fNeutron_Corrected_MomY_Col;
extern double fNeutron_Corrected_MomZ_Col;
extern double fNeutron_Delta_Mom_Col;
extern double fNeutron_Corrected_Energy_Col_GeV;
extern double fNeutron_Corrected_Mom_Col_GeV;
extern double fNeutron_Corrected_MomX_Col_GeV;
extern double fNeutron_Corrected_MomY_Col_GeV;
extern double fNeutron_Corrected_MomZ_Col_GeV;
extern double fNeutron_Delta_Mom_Col_GeV;

extern double fRecoilProton_Energy_RF;
extern double fRecoilProton_Mom_RF;
extern double fRecoilProton_MomX_RF;
extern double fRecoilProton_MomY_RF;
extern double fRecoilProton_MomZ_RF;
extern double fRecoilProton_Energy_RF_GeV;
extern double fRecoilProton_Mom_RF_GeV;
extern double fRecoilProton_MomX_RF_GeV;
extern double fRecoilProton_MomY_RF_GeV;
extern double fRecoilProton_MomZ_RF_GeV;
extern double fRecoilProton_Theta_RF;
extern double fRecoilProton_Phi_RF;

extern double fRecoilProton_Targ_Thickness;
extern double fRecoilProton_Targ_Thickness_RadLen;
extern double fRecoilProton_Targ_BT;
extern double fRecoilProton_Targ_Bremss_Loss;
extern double fRecoilProton_Targ_Ion_Loss;
extern double fRecoilProton_Air_Thickness;
extern double fRecoilProton_Air_Thickness_RadLen;
extern double fRecoilProton_Air_BT;
extern double fRecoilProton_Air_Bremss_Loss;
extern double fRecoilProton_Air_Ion_Loss;
extern double fRecoilProton_Theta_Col;
extern double fRecoilProton_Phi_Col;
extern double fRecoilProton_Energy_Col;
extern double fRecoilProton_Mom_Col;
extern double fRecoilProton_MomZ_Col;
extern double fRecoilProton_MomX_Col;
extern double fRecoilProton_MomY_Col;
extern double fRecoilProton_Energy_Col_GeV;
extern double fRecoilProton_Mom_Col_GeV;
extern double fRecoilProton_MomX_Col_GeV;
extern double fRecoilProton_MomY_Col_GeV;
extern double fRecoilProton_MomZ_Col_GeV;
extern double fRecoilProton_Corrected_Theta_Col;
extern double fRecoilProton_Corrected_Phi_Col;
extern double fRecoilProton_Corrected_Energy_Col;
extern double fRecoilProton_Corrected_Mom_Col;
extern double fRecoilProton_Corrected_MomX_Col;
extern double fRecoilProton_Corrected_MomY_Col;
extern double fRecoilProton_Corrected_MomZ_Col;
extern double fRecoilProton_Delta_Mom_Col;
extern double fRecoilProton_Corrected_Energy_Col_GeV;
extern double fRecoilProton_Corrected_Mom_Col_GeV;
extern double fRecoilProton_Corrected_MomX_Col_GeV;
extern double fRecoilProton_Corrected_MomY_Col_GeV;
extern double fRecoilProton_Corrected_MomZ_Col_GeV;
extern double fRecoilProton_Delta_Mom_Col_GeV;

extern double fSSAsym;
extern double fSineAsym;
extern double fInvariantDif;
extern double fT_GeV;
extern double fProton_Kin_Col;
extern double fQsq_Value;
extern double fQsq_Dif;
extern double fQsq_GeV;
extern double fQsq;
extern double fW_GeV_Col;
extern double fW_Col;
extern double fW;
extern double fW_GeV;
extern double fW_Prime_GeV;
extern double fW_Corrected_Prime_GeV;
extern double fWSq;
extern double fWSq_GeV;
extern double fWSq_PiN;
extern double fWSq_PiN_GeV;
extern double fWSq_Top_PiN_GeV;
extern double fWSq_Bot_PiN_GeV;

extern double fElec_ScatElec_Theta_RF;
extern double fScatElec_Cone_Phi_RF;
extern double fScatElec_Theta_RF;
extern double fScatElec_Phi_RF;
extern double fScatElec_Mom_RF;
extern double fScatElec_Energy_RF;
extern double fScatElec_MomX_RF;
extern double fScatElec_MomZ_RF;
extern double fScatElec_MomY_RF;
extern double fScatElec_Energy_RF_GeV;
extern double fScatElec_Mom_RF_GeV;
extern double fScatElec_MomX_RF_GeV;
extern double fScatElec_MomY_RF_GeV;
extern double fScatElec_MomZ_RF_GeV;

extern double fElectron_Theta_RF;
extern double fElectron_Phi_RF;
extern double fElectron_Energy_RF;
extern double fElectron_Mom_RF;
extern double fElectron_MomX_RF;
extern double fElectron_MomZ_RF;
extern double fElectron_MomY_RF;
extern double fElectron_Energy_RF_GeV;
extern double fElectron_Mom_RF_GeV;
extern double fElectron_MomX_RF_GeV;
extern double fElectron_MomZ_RF_GeV;
extern double fElectron_MomY_RF_GeV;

extern double fPhoton_Energy_RF_GeV;
extern double fPhoton_Mom_RF_GeV;
extern double fPhoton_Energy_RF;
extern double fPhoton_Mom_RF;

extern double fProton_Energy_CM;
extern double fProton_Mom_CM;
extern double fProton_Energy_CM_GeV;
extern double fProton_Mom_CM_GeV;
extern double fPhoton_Energy_CM;
extern double fPhoton_Mom_CM;
extern double fPhoton_Energy_CM_GeV;
extern double fPhoton_Mom_CM_GeV;
extern double fPion_Theta_CM;
extern double fPion_Phi_CM;
extern double fPion_Energy_CM;
extern double fPion_Mom_CM;
extern double fPion_Energy_CM_GeV;
extern double fPion_Mom_CM_GeV;
extern double fNeutron_Theta_CM;
extern double fNeutron_Phi_CM;
extern double fNeutron_Energy_CM;
extern double fNeutron_Energy_CM_GeV;
extern double fNeutron_Mom_CM;
extern double fNeutron_Mom_CM_GeV;

extern double fBeta_CM_RF;
extern double fGamma_CM_RF;

extern double fPhoton_MomZ_RF;
extern double fPhoton_MomX_RF;
extern double fPhoton_MomY_RF;
extern double fPhoton_Theta_RF;
extern double fPhoton_Phi_RF;
extern double fPion_Energy_RF;
extern double fPion_Energy_RF_GeV;
extern double fPiqVec_Theta_RF;
extern double fPion_Mom_RF;
extern double fPion_Mom_RF_GeV;
extern double fPion_MomX_RF;
extern double fPion_MomY_RF;
extern double fPion_MomZ_RF;
extern double fPion_Theta_RF;
extern double fPion_Phi_RF;
extern double fPion_MomX_RF_GeV;
extern double fPion_MomY_RF_GeV;
extern double fPion_MomZ_RF_GeV;

extern double fT_Para;
extern double fT_Para_GeV;
extern double fT;
extern double fEpsilon;
extern double fx;
extern double fy;
extern double fz;
extern double fNeutron_Energy_RF;
extern double fNeutron_Energy_RF_GeV;
extern double fNeutron_Mom_RF;
extern double fNeutron_Mom_RF_GeV;
extern double fNeutron_qVec_Theta_RF;
extern double fNeutron_MomX_RF;
extern double fNeutron_MomY_RF;
extern double fNeutron_MomZ_RF;
extern double fNeutron_Theta_RF;
extern double fNeutron_Phi_RF;
extern double fPhoton_MomX_RF_GeV;
extern double fPhoton_MomY_RF_GeV;
extern double fPhoton_MomZ_RF_GeV;
extern double fNeutron_MomX_RF_GeV;
extern double fNeutron_MomY_RF_GeV;
extern double fNeutron_MomZ_RF_GeV;

extern double fPhoton_Theta_Col;
extern double fPhoton_Phi_Col;
extern double fPhoton_Energy_Col;
extern double fPhoton_Mom_Col;
extern double fPhoton_MomX_Col;
extern double fPhoton_MomZ_Col;
extern double fPhoton_MomY_Col;
extern double fPhoton_Energy_Col_GeV;
extern double fPhoton_Mom_Col_GeV;
extern double fPhoton_MomX_Col_GeV;
extern double fPhoton_MomZ_Col_GeV;
extern double fPhoton_MomY_Col_GeV;

extern double fPhoton_Corrected_Theta_Col;
extern double fPhoton_Corrected_Phi_Col;
extern double fPhoton_Corrected_Energy_Col;
extern double fPhoton_Corrected_Mom_Col;
extern double fPhoton_Corrected_MomX_Col;
extern double fPhoton_Corrected_MomZ_Col;
extern double fPhoton_Corrected_MomY_Col;
extern double fPhoton_Corrected_Energy_Col_GeV;
extern double fPhoton_Corrected_Mom_Col_GeV;
extern double fPhoton_Corrected_MomX_Col_GeV;
extern double fPhoton_Corrected_MomZ_Col_GeV;
extern double fPhoton_Corrected_MomY_Col_GeV;

extern double fQsq_Corrected_GeV;
extern double fQsq_Corrected;
extern double fW_Corrected;
extern double fW_Corrected_GeV;
extern double fT_Corrected;
extern double fT_Corrected_GeV;
extern double fx_Corrected;
extern double fy_Corrected;
extern double fz_Corrected;

extern double fWFactor;
extern double fA;
extern double fFlux_Factor_Col;
extern double fFlux_Factor_RF;
extern double fJacobian_CM;
extern double fJacobian_CM_RF;
extern double fJacobian_CM_Col;
extern double fZASig_T;
extern double fZASig_L;
extern double fZASig_LT;
extern double fZASig_TT;
extern double ftestsig;
extern double fZASig_L2;

extern double fZASigma_UU;
extern double fRorySigma_UT;
extern double fSigma_Col;
extern double fSigma_UUPara;
extern double fSig_VR;
extern double fSig_L;
extern double fSig_T;

extern double fSig_fpi_6GeV;

extern double fSigmaPhiS;
extern double fSigmaPhi_Minus_PhiS;
extern double fSigma2Phi_Minus_PhiS;
extern double fSigma3Phi_Minus_PhiS;
extern double fSigmaPhi_Plus_PhiS;
extern double fSigma2Phi_Plus_PhiS;
extern double fSig_Phi_Minus_PhiS;
extern double fSig_PhiS;
extern double fSig_2Phi_Minus_PhiS;
extern double fSig_Phi_Plus_PhiS;
extern double fSig_3Phi_Minus_PhiS;
extern double fSig_2Phi_Plus_PhiS;
extern double fEventWeight;
extern double fEventWeightMax;
extern double fEventWeightCeil;  // SJDK 11/05/21 - This is the maximum value found with the old method that is used to get the new unit weight
extern double fEventWeightRn;  // SJDK 11/05/21 - Random number to compare determined weight to
extern double fZAWFactor;
extern double fRR;
extern double fPhaseSpaceWeight;
extern double fPhaseShiftWeight;
extern double fWilliamsWeight;
extern double fDedrickWeight;
extern double fCatchenWeight;
extern double fPhi;
extern double fPhiS;
extern double fPhi_Corrected;
extern double fPhiS_Corrected;

extern double fElectron_Mom_Sq_RF;
extern double fElectron_Mom_Sq_Col;
extern double fProton_Mom_Sq_Col;
extern double fProton_Mom_Sq_CM;
extern double fProton_Mom_Sq_RF;
extern double fPhoton_Mom_Sq_Col;
extern double fPhoton_Mom_Sq_CM;
extern double fPhoton_Mom_Sq_RF;
extern double fPion_Mom_Sq_Col;
extern double fPion_Mom_Sq_CM;
extern double fPion_Mom_Sq_RF;
extern double fNeutron_Mom_Sq_Col;
extern double fNeutron_Mom_Sq_CM;
extern double fNeutron_Mom_Sq_RF;
extern double fScatElec_Mom_Sq_Col;
extern double fScatElec_Mom_Sq_RF;

extern double fAsymPhiMinusPhi_S;
extern double fAsymPhi_S;
extern double fAsym2PhiMinusPhi_S;
extern double fAsymPhiPlusPhi_S;
extern double fAsym3PhiMinusPhi_S;
extern double fAsym2PhiPlusPhi_S;

extern double fTerm_PhiMinusPhi_S;
extern double fTerm_Phi_S;
extern double fTerm_2PhiMinusPhi_S;
extern double fTerm_PhiPlusPhi_S;
extern double fTerm_3PhiMinusPhi_S;
extern double fTerm_2PhiPlusPhi_S;

extern double fAsymPhiMinusPhi_S_Col;
extern double fAsymPhi_S_Col;
extern double fAsym2PhiMinusPhi_S_Col;
extern double fAsymPhiPlusPhi_S_Col;
extern double fAsym3PhiMinusPhi_S_Col;
extern double fAsym2PhiPlusPhi_S_Col;

extern double fTerm_PhiMinusPhi_S_Col;
extern double fTerm_Phi_S_Col;
extern double fTerm_2PhiMinusPhi_S_Col;
extern double fTerm_PhiPlusPhi_S_Col;
extern double fTerm_3PhiMinusPhi_S_Col;
extern double fTerm_2PhiPlusPhi_S_Col;

extern double fPhi_Pion_LeptonPlane_RF;
extern double fCos_Phi_Pion_LeptonPlane_RF;
extern double fSin_Phi_Pion_LeptonPlane_RF;
extern double fPhi_TargPol_LeptonPlane_RF;
extern double fCos_Phi_TargPol_LeptonPlane_RF;
extern double fSin_Phi_TargPol_LeptonPlane_RF;
extern double fTheta_Pion_Photon_RF;
extern double fPhi_Pion_LeptonPlane_Col;
extern double fCos_Phi_Pion_LeptonPlane_Col;
extern double fSin_Phi_Pion_LeptonPlane_Col;
extern double fPhi_TargPol_LeptonPlane_Col;
extern double fCos_Phi_TargPol_LeptonPlane_Col;
extern double fSin_Phi_TargPol_LeptonPlane_Col;
extern double fTheta_Pion_Photon_Col;

extern double fZASigma_UU_Col;
extern double fRorySigma_UT_Col;
extern double fSig_Phi_Minus_PhiS_Col;
extern double fSig_PhiS_Col;
extern double fSig_2Phi_Minus_PhiS_Col;
extern double fSig_Phi_Plus_PhiS_Col;
extern double fSig_3Phi_Minus_PhiS_Col;
extern double fSig_2Phi_Plus_PhiS_Col;

extern double fepi1;
extern double fepi2;
extern double fradical;

extern double fMomentum[300];
extern double fProb[300];

extern double conserve;      // 16/06/21 AU -> New Variables for conservation law checks. This is the number that PASS both conservation check
extern double ene; // The number that FAIL due to failing energy conservation check ONLY
extern double mom; // The number that FAIL due to failing the momentum conservation check ONLY
extern double ene_mom; // The number that FAIL BOTH energy and momentum check

extern double mom_px; // Number that fail due to px only
extern double mom_py; // Number that fail due to py only
extern double mom_pz; // Number that fail due to pz only
extern double mom_pxpy; // Number that fail due to px and py
extern double mom_pxpz; // Number that fail due to px and pz
extern double mom_pypz; // Number that fail due to py and pz
extern double mom_pxpypz; // Number that fail due to px, py and pz

// 27/01/22 - Love Preet - Adding in vector of cross section parameters
extern vector<vector<vector<vector<double>>>> SigPar;

// Love Preet - Added for phase space factor calculations
extern int psf_steps;
extern double psf_ScatElec_E_Stepsize;
extern double psf_ScatElec_Theta_Stepsize;
extern double psf_ScatElec_Phi_Stepsize;
extern double psf_ScatElec_E;
extern double psf_ScatElec_Theta;
extern double psf_ScatElec_Phi;
extern double psf_ScalElec_Mom;
extern double psf_Ejec_Theta_Stepsize;
extern double psf_Ejectile_Theta;
extern double psf_Ejectile_Phi;
extern double psf_Q2, psf_W, psf_W2, psf_t;
extern double psf_ScatElec_Theta_max, psf_ScatElec_Theta_min;
extern double psf_ScatElec_E_max, psf_ScatElec_E_min;
extern double psf_Ejectile_Theta_max, psf_Ejectile_Theta_min;

// Love Preet - Added for actual phase space factor calculations
extern double fScatElec_Energy_Col_max;
extern double fScatElec_Energy_Col_min; 
extern double fScatElec_Theta_Col_max; 
extern double fScatElec_Theta_Col_min;
extern double f_Ejectile_Theta_Col_max;
extern double f_Ejectile_Theta_Col_min;
extern double fPSF_org;

//Love Preet - Added to be stored in the root tree
extern double scat_e_px;
extern double scat_e_py;
extern double scat_e_pz;
extern double scat_e_E;
extern double ejec_px;
extern double ejec_py;
extern double ejec_pz;
extern double ejec_E;
extern double rclH_px;
extern double rclH_py;
extern double rclH_pz;
extern double rclH_E;

#endif
