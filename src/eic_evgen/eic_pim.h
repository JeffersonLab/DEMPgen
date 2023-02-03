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
  	int CheckLaws(TLorentzVector P_E0, TLorentzVector P_t, TLorentzVector P_e, TLorentzVector P_pim, TLorentzVector P_pro, double fdiff_E);
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
extern std::string gDet_location;
extern std::string gOutputType;
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

extern unsigned long long int fNaN;
extern unsigned long long int fConserve;

extern unsigned long long int fNWeightUnphys;
extern unsigned long long int fNWeightReject;

extern unsigned long long int fLundRecorded;
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
extern double fScatElec_Theta_I;
extern double fScatElec_Theta_F;
extern double fPion_Theta_I; // SJDK 19/12/22 - These should be removed in future, specific to pion reaction cases. Should be generic MesonX
extern double fPion_Theta_F;
extern double fEjectileX_Theta_I;
extern double fEjectileX_Theta_F;
extern double fScatElec_E_Hi;
extern double fScatElec_E_Lo;
extern double fPSF;

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

extern double conserve;      // 16/06/21 AU -> New Variables for conservation law checks
extern double ene;
extern double mom;

// 27/01/22 - Love Preet - Adding in vector of cross section parameters
extern vector<vector<vector<vector<double>>>> SigPar;

//extern double fProb[300] = {    
//6.03456,    6.02429,    6.01155,    5.99636,    5.97873,    5.95869,    5.93626,    5.91147,    5.88435,    5.85493,
//5.82325,    5.78935,    5.75326,    5.71504,    5.67472,    5.63235,    5.58799,    5.54169,     5.4935,    5.44347,		   
//5.39167,    5.33816,    5.28299,    5.22623,    5.16794,    5.10818,    5.04703,    4.98455,    4.92081,    4.85588,
//4.78982,    4.71692,    4.63621,    4.55583,    4.47582,    4.39621,    4.31702,    4.23828,    4.16002,    4.08227,
//4.00506,     3.9284,    3.85233,    3.77686,    3.70202,    3.62783,    3.55432,    3.48149,    3.40937,    3.33798,
//3.26733,    3.19745,    3.12834,    3.06002,    2.99251,    2.92581,    2.85995,    2.79493,    2.73075,    2.66744,
//2.605,      2.54344,    2.48276,    2.41728,    2.35244,    2.28922,    2.22759,    2.16751,    2.10895,    2.05186,
//1.99621,    1.94198,    1.88913,    1.83762,    1.78743,    1.73851,    1.69086,    1.64442,    1.59918,    1.55511,			   
//1.51217,    1.47035,    1.42961,    1.38993,    1.35128,    1.31364,    1.27698,    1.24129,    1.20653,    1.17269,			   
//1.13973,    1.10765,    1.07642,    1.04601,    1.01626,   0.986934,   0.958443,   0.930759,   0.903861,   0.877727,			   
//0.852335,   0.827665,   0.803696,   0.780409,   0.757785,   0.735806,   0.714452,   0.693708,   0.673555,   0.653978,			   
//0.63496,   0.616485,   0.598538,   0.581105,    0.56417,   0.547721,   0.531743,   0.516223,   0.501148,   0.486506,			   
//0.472284,    0.45847,   0.445054,   0.432024,   0.419368,   0.407077,   0.395125,   0.383513,   0.372237,   0.361289,			   
//0.350658,   0.340336,   0.330314,   0.320583,   0.311135,   0.301962,   0.293056,   0.284409,   0.276014,   0.267863,			   
//0.25995,   0.252268,   0.244809,   0.237569,   0.230539,   0.223715,   0.217091,    0.21066,   0.204417,   0.198357,			   
//0.192474,   0.186763,   0.181219,   0.175838,   0.170615,   0.165545,   0.160623,   0.155843,   0.151185,   0.146666,			   
//0.142282,   0.138028,   0.133902,   0.129898,   0.126013,   0.122245,   0.118588,   0.115041,   0.111599,   0.108261,			   
//0.105021,   0.101879,  0.0988298,  0.0958719,  0.0930022,  0.0902182,  0.0875173,   0.084897,  0.0823549,  0.0798886,			   
//0.0774961,  0.0751749,  0.0729231,  0.0707385,  0.0686192,  0.0665631,  0.0645685,  0.0626335,  0.0607563,  0.0589545,			   
//0.057219,  0.0555322,   0.053893,  0.0523001,  0.0507522,  0.0492481,  0.0477867,  0.0463667,  0.0449872,  0.0436468,			   
//0.0423448,  0.0410798,  0.0398511,  0.0386576,  0.0374982,  0.0363722,  0.0352785,  0.0342164,  0.0331849,  0.0321831,			   
//0.0312104,  0.0302658,  0.0293486,  0.0284581,  0.0275935,   0.026754,  0.0259391,  0.0251479,  0.0243799,  0.0236344,			   
//0.0229107,  0.0221901,  0.0214923,  0.0208167,  0.0201626,  0.0195293,  0.0189162,  0.0183226,  0.0177478,  0.0171913,			   
//0.0166525,  0.0161308,  0.0156257,  0.0151366,  0.0146629,  0.0142044,  0.0137603,  0.0133303,  0.0129139,  0.0125107,			   
//0.0121203,  0.0117421,   0.011376,  0.0110214,   0.010678,  0.0103455,  0.0100234, 0.00971153, 0.00940947, 0.00911693,			   
//0.00883361,  0.0085592, 0.00829385, 0.00803735, 0.00778883, 0.00754804, 0.00731474,  0.0070887, 0.00686967, 0.00665746,			   
//0.00645184, 0.00625261, 0.00605957, 0.00587252, 0.00569128, 0.00551567, 0.00534551, 0.00518063, 0.00502086, 0.00486605,			   
//0.00471604, 0.00457069, 0.00442984, 0.00429336,  0.0041611, 0.00403295, 0.00390877, 0.00378843, 0.00367182, 0.00355882,			   
//0.00344932, 0.00334321, 0.00324038, 0.00314073, 0.00304505,  0.0029524, 0.00286252, 0.00277533, 0.00269076, 0.00260872,			   
//0.00252913, 0.00245194, 0.00237706, 0.00230444, 0.00223399, 0.00216566, 0.00209939, 0.00203512, 0.00197277, 0.00191231 };
//

#endif
