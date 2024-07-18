///*--------------------------------------------------*/
/// eic_pim.cc:
/// Originally did not exist
/// Date: 2015-2018
///
///*--------------------------------------------------*/
/// Modifier: Wenliang (Bill) Li
/// Date: Feb 24 2020
/// Email: wenliang.billlee@gmail.com
///
/// Comment: Feb 24, 2020: all declearation of the global variables 
//            and all functions for pim class are moved there 

#include "eic_pim.h"

using std::vector;
using namespace std;

//TRandom2 *fRandom;                    
TRandom3 *fRandom;                    

TFile *f;                    
TTree *t1;

int gKinematics_type;
bool gPi0_decay;
bool UseSolve;
bool gROOTOut; // SJDK 30/04/24 - Added boolean to enable/disable ROOTfile output
string gDet_location;
string gOutputType; // SJDK 12/01/22 - Added output type as a variable you can specify in the .json file
string gBeamPart; // SJDK 12/01/22 - Added output type as a variable you can specify in the .json file
float fProton_incidence_phi;

TString gfile_name;

int fSeed;

bool allset, print, kCalcFermi, kCalcBremss, kCalcIon, kCalcBremssEle, kCalcIonEle, kSConserve, kFSI, kMSele, kMS;
int fWLessShell, fWLess1P9, fSDiff;

//long int fNEvents, fNRecorded, fNGenerated, fWSqNeg, fNSigmaNeg, fNWeightUnphys, fNWeightReject, fLundRecorded, fNFile; 

unsigned long long int fNEvents, fNRecorded, fNGenerated, fWSqNeg, fNSigmaNeg, fNaN, fConserve, fNWeightUnphys, fNWeightReject, fLundRecorded, fNFile, fSolveEvents_0Sol, fSolveEvents_1Sol, fSolveEvents_2Sol, fNWeightNeg; // Love Preet - Added the fNWeightNeg to get number of negative weights

// SJDK 03/04/23 - Added in Qsq Min/Max, W Min/Max and t max (06/09/23)
// 13/09/23 - SJDK - New generic HBeam value (rather than proton beam)
double fK, fm, fElectron_Kin_Col_GeV, fElectron_Kin_Col, fRand, fLumi, fuBcm2, fPI, fDEG2RAD, fRAD2DEG, fEBeam, fPBeam, fHBeam, fScatElec_Theta_I, fScatElec_Theta_F, fPion_Theta_I, fPion_Theta_F, fEjectileX_Theta_I, fEjectileX_Theta_F, fScatElec_E_Hi, fScatElec_E_Lo, fPSF, fQsq_Min, fQsq_Max, fW_Min, fW_Max, fT_Max; 

double fOmega_Theta_I, fOmega_Theta_F, fOmega_Theta_Col, fOmega_Phi_Col;

double fDiff_E, conserve, ene, mom, ene_mom, mom_px, mom_py, mom_pz, mom_pxpy, mom_pxpz, mom_pypz, mom_pxpypz; // 18/06/21 AU -> New variables to count envents passing/not passing conservation laws

double fMandSConserve, fTop_Pion_Mom, fBot_Pion_Mom, fPion_Mom_Same, fEnergyConserve, fXMomConserve, fYMomConserve, fZMomConserve, fXMomConserve_RF, fYMomConserve_RF, fZMomConserve_RF, fEnergyConserve_RF; 

double fDiff, fRatio, fPion_Alpha, fPion_Beta, fS_I_RF, fS_F_RF, fS_I_Col, fS_F_Col, fS_I_RF_GeV, fS_F_RF_GeV, fS_I_Col_GeV, fS_F_Col_GeV; 

double fProton_Energy_Col, fProton_Mom_Col, fProton_Theta_Col, fProton_Phi_Col, fProton_MomZ_Col, fProton_MomX_Col, fProton_MomY_Col, fProton_Energy_Col_GeV, fProton_Mom_Col_GeV, fProton_MomX_Col_GeV, fProton_MomY_Col_GeV, fProton_MomZ_Col_GeV; 

double fFSIProton_Energy_Col, fFSIProton_Mom_Col, fFSIProton_Theta_Col, fFSIProton_Phi_Col, fFSIProton_MomZ_Col, fFSIProton_MomX_Col, fFSIProton_MomY_Col, fFSIProton_Energy_Col_GeV, fFSIProton_Mom_Col_GeV, fFSIProton_MomX_Col_GeV, fFSIProton_MomY_Col_GeV, fFSIProton_MomZ_Col_GeV;

double fTarget_Energy_Col, fTarget_Mom_Col, fTarget_Theta_Col, fTarget_Phi_Col, fTarget_MomZ_Col, fTarget_MomX_Col, fTarget_MomY_Col, fTarget_Energy_Col_GeV, fTarget_Mom_Col_GeV, fTarget_MomX_Col_GeV, fTarget_MomY_Col_GeV, fTarget_MomZ_Col_GeV;

double fTarget_Pol0_Col, fTarget_PolX_Col, fTarget_PolY_Col, fTarget_PolZ_Col, fTarget_Pol0_RF, fTarget_PolX_RF, fTarget_PolY_RF, fTarget_PolZ_RF;

double fBetaX_Col_RF, fBetaY_Col_RF, fBetaZ_Col_RF, fBeta_Col_RF, fGamma_Col_RF;

double fProton_MomX_RF, fProton_MomY_RF, fProton_MomZ_RF, fProton_Mom_RF, fProton_Energy_RF, fProton_MomX_RF_GeV, fProton_MomY_RF_GeV, fProton_MomZ_RF_GeV, fProton_Mom_RF_GeV, fProton_Energy_RF_GeV;

double fScatElec_Angle, fScatElec_Alpha_RF, fScatElec_Beta_RF;

double fVertex_X, fVertex_Y, fVertex_Z, fProton_Kin_Col_GeV, fElectron_Mass, fElectron_Mass_GeV, fProton_Mass, fProton_Mass_GeV, fNeutron_Mass, fNeutron_Mass_GeV, fPion_Mass, fPion_Mass_GeV, fPiion_Phi, fAlpha, fPi, fMom_Ratio, fMom_Dif, fPionEnergyCMLess, fSNotEqual, fMode_Epsi, fRecoilProton_Mass, fRecoilProton_Mass_GeV;

double fOmega_Mass, fOmega_Mass_GeV; 

double f_Scat_hadron_Mass, f_Scat_hadron_Mass_GeV;

double fKaon_Mass, fKaon_Mass_GeV, fLambda_Mass, fLambda_Mass_GeV, fSigma_Mass, fSigma_Mass_GeV;

double fElectron_Energy_Col, fElectron_MomZ_Col, fElectron_MomX_Col, fElectron_MomY_Col, fElectron_Theta_Col, fElectron_Phi_Col, fElectron_Mom_Col;

double fElectron_MS_Energy_Col, fElectron_MS_MomZ_Col, fElectron_MS_MomX_Col, fElectron_MS_MomY_Col, fElectron_MS_Theta_Col, fElectron_MS_Phi_Col, fElectron_MS_Mom_Col;

double fElectron_Energy_Col_GeV, fElectron_Mom_Col_GeV, fElectron_MomX_Col_GeV, fElectron_MomY_Col_GeV, fElectron_MomZ_Col_GeV, fElectronEnergyLess, fElectronThetaLess, fRadiation_Lenght_Air;

double fElectron_Targ_Thickness, fElectron_Targ_Thickness_RadLen, fElectron_Targ_BT, fElectron_Targ_Bremss_Loss, fElectron_Targ_Ion_Loss, fElectron_TargWindow_Bremss_Loss, fElectron_TargWindow_Ion_Loss;

double fElectron_Air_Thickness, fElectron_Air_Thickness_RadLen, fElectron_Air_BT, fElectron_Air_Bremss_Loss, fElectron_Air_Ion_Loss, fElectron_Corrected_Energy_Col, fElectron_Corrected_Mom_Col, fElectron_Corrected_MomX_Col, fElectron_Corrected_MomY_Col, fElectron_Corrected_MomZ_Col, fElectron_Corrected_Theta_Col, fElectron_Corrected_Phi_Col, fElectron_Delta_Mom_Col, fElectron_Corrected_Energy_Col_GeV, fElectron_Corrected_Mom_Col_GeV, fElectron_Corrected_MomX_Col_GeV, fElectron_Corrected_MomY_Col_GeV, fElectron_Corrected_MomZ_Col_GeV, fElectron_Delta_Mom_Col_GeV;

double fScatElec_MS_Energy_Col, fScatElec_MS_MomZ_Col, fScatElec_MS_MomX_Col, fScatElec_MS_MomY_Col, fScatElec_MS_Theta_Col, fScatElec_MS_Phi_Col, fScatElec_MS_Mom_Col;

double fScatElec_Energy_Col, fScatElec_MomZ_Col, fScatElec_MomX_Col, fScatElec_MomY_Col, fScatElec_Theta_Col, fScatElec_Phi_Col, fScatElec_Mom_Col, fScatElec_Energy_Col_GeV, fScatElec_Mom_Col_GeV, fScatElec_MomX_Col_GeV, fScatElec_MomY_Col_GeV, fScatElec_MomZ_Col_GeV, fScatElecEnergyLess, fScatElecThetaLess, fScatElec_Targ_Thickness, fScatElec_Targ_Thickness_RadLen, fScatElec_Targ_BT, fScatElec_Targ_Bremss_Loss, fScatElec_Targ_Ion_Loss, fScatElec_Air_Thickness, fScatElec_Air_Thickness_RadLen, fScatElec_Air_BT, fScatElec_Air_Bremss_Loss, fScatElec_Air_Ion_Loss, fScatElec_Corrected_Energy_Col, fScatElec_Corrected_Mom_Col, fScatElec_Corrected_MomX_Col, fScatElec_Corrected_MomY_Col, fScatElec_Corrected_MomZ_Col, fScatElec_Corrected_Theta_Col, fScatElec_Corrected_Phi_Col, fScatElec_Delta_Mom_Col, fScatElec_Corrected_Energy_Col_GeV, fScatElec_Corrected_Mom_Col_GeV, fScatElec_Corrected_MomX_Col_GeV, fScatElec_Corrected_MomY_Col_GeV, fScatElec_Corrected_MomZ_Col_GeV, fScatElec_Delta_Mom_Col_GeV, fScatElec_TargWindow_Bremss_Loss, fScatElec_TargWindow_Ion_Loss, fTargWindow_Thickness, fTargWindow_Thickness_RadLen, fTargWindow_BT;

double fPion_TargWindow_Ion_Loss, fPion_Targ_Thickness, fPion_Targ_Thickness_RadLen, fPion_Targ_BT, fPion_Targ_Bremss_Loss, fPion_Targ_Ion_Loss, fPion_Air_Thickness, fPion_Air_Thickness_RadLen, fPion_Air_BT, fPion_Air_Bremss_Loss, fPion_Air_Ion_Loss;

double fPion_MS_Energy_Col, fPion_MS_MomZ_Col, fPion_MS_MomX_Col, fPion_MS_MomY_Col, fPion_MS_Theta_Col, fPion_MS_Phi_Col, fPion_MS_Mom_Col;

double fPion_Theta_Col, fPion_Phi_Col, fPion_Energy_Col, fPion_Mom_Col, fPion_MomZ_Col, fPion_MomX_Col, fPion_MomY_Col, fPion_Energy_Col_GeV, fPion_Mom_Col_GeV, fPion_MomX_Col_GeV, fPion_MomY_Col_GeV, fPion_MomZ_Col_GeV;

double fKaon_Theta_Col, fKaon_Phi_Col, fKaon_Energy_Col, fKaon_Mom_Col, fKaon_MomZ_Col, fKaon_MomX_Col, fKaon_MomY_Col, fKaon_Energy_Col_GeV, fKaon_Mom_Col_GeV, fKaon_MomX_Col_GeV, fKaon_MomY_Col_GeV, fKaon_MomZ_Col_GeV;

double fScathad_Theta_Col, fScathad_Phi_Col, fScathad_Energy_Col, fScathad_Mom_Col, fScathad_MomZ_Col, fScathad_MomX_Col, fScathad_MomY_Col, fScathad_Energy_Col_GeV, fScathad_Mom_Col_GeV, fScathad_MomX_Col_GeV, fScathad_MomY_Col_GeV, fScathad_MomZ_Col_GeV;

double fPion_FSI_Theta_Col, fPion_FSI_Phi_Col, fPion_FSI_Energy_Col, fPion_FSI_Mom_Col, fPion_FSI_MomZ_Col, fPion_FSI_MomX_Col, fPion_FSI_MomY_Col, fPion_FSI_Energy_Col_GeV, fPion_FSI_Mom_Col_GeV, fPion_FSI_MomX_Col_GeV, fPion_FSI_MomY_Col_GeV, fPion_FSI_MomZ_Col_GeV;

double fPion_Corrected_Theta_Col, fPion_Corrected_Phi_Col, fPion_Corrected_Energy_Col, fPion_Corrected_Mom_Col, fPion_Corrected_MomX_Col, fPion_Corrected_MomY_Col, fPion_Corrected_MomZ_Col, fPion_Delta_Mom_Col, fPion_Corrected_Energy_Col_GeV, fPion_Corrected_Mom_Col_GeV, fPion_Corrected_MomX_Col_GeV, fPion_Corrected_MomY_Col_GeV, fPion_Corrected_MomZ_Col_GeV, fPion_Delta_Mom_Col_GeV;

double fNeutron_MS_Energy_Col, fNeutron_MS_MomZ_Col, fNeutron_MS_MomX_Col, fNeutron_MS_MomY_Col, fNeutron_MS_Theta_Col, fNeutron_MS_Phi_Col, fNeutron_MS_Mom_Col;

double fNeutron_TargWindow_Ion_Loss, fNeutron_Targ_Thickness, fNeutron_Targ_Thickness_RadLen, fNeutron_Targ_BT, fNeutron_Targ_Bremss_Loss, fNeutron_Targ_Ion_Loss, fNeutron_Air_Thickness, fNeutron_Air_Thickness_RadLen, fNeutron_Air_BT, fNeutron_Air_Bremss_Loss, fNeutron_Air_Ion_Loss, fNeutron_Theta_Col, fNeutron_Phi_Col, fNeutron_Energy_Col, fNeutron_Mom_Col, fNeutron_MomZ_Col, fNeutron_MomX_Col, fNeutron_MomY_Col, fNeutron_Energy_Col_GeV, fNeutron_Mom_Col_GeV, fNeutron_MomX_Col_GeV, fNeutron_MomY_Col_GeV, fNeutron_MomZ_Col_GeV, fNeutron_Corrected_Theta_Col, fNeutron_Corrected_Phi_Col, fNeutron_Corrected_Energy_Col, fNeutron_Corrected_Mom_Col, fNeutron_Corrected_MomX_Col, fNeutron_Corrected_MomY_Col, fNeutron_Corrected_MomZ_Col, fNeutron_Delta_Mom_Col, fNeutron_Corrected_Energy_Col_GeV, fNeutron_Corrected_Mom_Col_GeV, fNeutron_Corrected_MomX_Col_GeV, fNeutron_Corrected_MomY_Col_GeV, fNeutron_Corrected_MomZ_Col_GeV, fNeutron_Delta_Mom_Col_GeV;

double fRecoilProton_Energy_RF, fRecoilProton_Mom_RF, fRecoilProton_MomX_RF, fRecoilProton_MomY_RF, fRecoilProton_MomZ_RF, fRecoilProton_Energy_RF_GeV, fRecoilProton_Mom_RF_GeV, fRecoilProton_MomX_RF_GeV, fRecoilProton_MomY_RF_GeV, fRecoilProton_MomZ_RF_GeV, fRecoilProton_Theta_RF, fRecoilProton_Phi_RF;

double fRecoilProton_Targ_Thickness, fRecoilProton_Targ_Thickness_RadLen, fRecoilProton_Targ_BT, fRecoilProton_Targ_Bremss_Loss, fRecoilProton_Targ_Ion_Loss, fRecoilProton_Air_Thickness, fRecoilProton_Air_Thickness_RadLen, fRecoilProton_Air_BT, fRecoilProton_Air_Bremss_Loss, fRecoilProton_Air_Ion_Loss, fRecoilProton_Theta_Col, fRecoilProton_Phi_Col, fRecoilProton_Energy_Col, fRecoilProton_Mom_Col, fRecoilProton_MomZ_Col, fRecoilProton_MomX_Col, fRecoilProton_MomY_Col, fRecoilProton_Energy_Col_GeV, fRecoilProton_Mom_Col_GeV, fRecoilProton_MomX_Col_GeV, fRecoilProton_MomY_Col_GeV, fRecoilProton_MomZ_Col_GeV, fRecoilProton_Corrected_Theta_Col, fRecoilProton_Corrected_Phi_Col, fRecoilProton_Corrected_Energy_Col, fRecoilProton_Corrected_Mom_Col, fRecoilProton_Corrected_MomX_Col, fRecoilProton_Corrected_MomY_Col, fRecoilProton_Corrected_MomZ_Col, fRecoilProton_Delta_Mom_Col, fRecoilProton_Corrected_Energy_Col_GeV, fRecoilProton_Corrected_Mom_Col_GeV, fRecoilProton_Corrected_MomX_Col_GeV, fRecoilProton_Corrected_MomY_Col_GeV, fRecoilProton_Corrected_MomZ_Col_GeV, fRecoilProton_Delta_Mom_Col_GeV;

double fSSAsym, fSineAsym, fInvariantDif, fT_GeV, fProton_Kin_Col, fQsq_Value, fQsq_Dif, fQsq_GeV, fQsq, fW_GeV_Col, fW_Col, fW, fW_GeV, fW_Prime_GeV, fW_Corrected_Prime_GeV, fWSq, fWSq_GeV, fWSq_PiN, fWSq_PiN_GeV, fWSq_Top_PiN_GeV, fWSq_Bot_PiN_GeV;

double fElec_ScatElec_Theta_RF, fScatElec_Cone_Phi_RF, fScatElec_Theta_RF, fScatElec_Phi_RF, fScatElec_Mom_RF, fScatElec_Energy_RF, fScatElec_MomX_RF, fScatElec_MomZ_RF, fScatElec_MomY_RF, fScatElec_Energy_RF_GeV, fScatElec_Mom_RF_GeV, fScatElec_MomX_RF_GeV, fScatElec_MomY_RF_GeV, fScatElec_MomZ_RF_GeV;

double fElectron_Theta_RF, fElectron_Phi_RF, fElectron_Energy_RF, fElectron_Mom_RF, fElectron_MomX_RF, fElectron_MomZ_RF, fElectron_MomY_RF, fElectron_Energy_RF_GeV, fElectron_Mom_RF_GeV, fElectron_MomX_RF_GeV, fElectron_MomZ_RF_GeV, fElectron_MomY_RF_GeV;

double fPhoton_Energy_RF_GeV, fPhoton_Mom_RF_GeV, fPhoton_Energy_RF, fPhoton_Mom_RF;

double fProton_Energy_CM, fProton_Mom_CM, fProton_Energy_CM_GeV, fProton_Mom_CM_GeV, fPhoton_Energy_CM, fPhoton_Mom_CM, fPhoton_Energy_CM_GeV, fPhoton_Mom_CM_GeV, fPion_Theta_CM, fPion_Phi_CM, fPion_Energy_CM, fPion_Mom_CM, fPion_Energy_CM_GeV, fPion_Mom_CM_GeV, fNeutron_Theta_CM, fNeutron_Phi_CM, fNeutron_Energy_CM, fNeutron_Energy_CM_GeV, fNeutron_Mom_CM, fNeutron_Mom_CM_GeV;

double fBeta_CM_RF, fGamma_CM_RF;

double fPhoton_MomZ_RF, fPhoton_MomX_RF, fPhoton_MomY_RF, fPhoton_Theta_RF, fPhoton_Phi_RF, fPion_Energy_RF, fPion_Energy_RF_GeV, fPiqVec_Theta_RF, fPion_Mom_RF, fPion_Mom_RF_GeV, fPion_MomX_RF, fPion_MomY_RF, fPion_MomZ_RF, fPion_Theta_RF, fPion_Phi_RF, fPion_MomX_RF_GeV, fPion_MomY_RF_GeV, fPion_MomZ_RF_GeV;

double fT_Para, fT_Para_GeV, fT, fEpsilon, fx, fy, fz, fNeutron_Energy_RF, fNeutron_Energy_RF_GeV, fNeutron_Mom_RF, fNeutron_Mom_RF_GeV, fNeutron_qVec_Theta_RF, fNeutron_MomX_RF, fNeutron_MomY_RF, fNeutron_MomZ_RF, fNeutron_Theta_RF, fNeutron_Phi_RF, fPhoton_MomX_RF_GeV, fPhoton_MomY_RF_GeV, fPhoton_MomZ_RF_GeV, fNeutron_MomX_RF_GeV, fNeutron_MomY_RF_GeV, fNeutron_MomZ_RF_GeV;

double fPhoton_Theta_Col, fPhoton_Phi_Col, fPhoton_Energy_Col, fPhoton_Mom_Col, fPhoton_MomX_Col, fPhoton_MomZ_Col, fPhoton_MomY_Col, fPhoton_Energy_Col_GeV, fPhoton_Mom_Col_GeV, fPhoton_MomX_Col_GeV, fPhoton_MomZ_Col_GeV, fPhoton_MomY_Col_GeV;

double fPhoton_Corrected_Theta_Col, fPhoton_Corrected_Phi_Col, fPhoton_Corrected_Energy_Col, fPhoton_Corrected_Mom_Col, fPhoton_Corrected_MomX_Col, fPhoton_Corrected_MomZ_Col, fPhoton_Corrected_MomY_Col, fPhoton_Corrected_Energy_Col_GeV, fPhoton_Corrected_Mom_Col_GeV, fPhoton_Corrected_MomX_Col_GeV, fPhoton_Corrected_MomZ_Col_GeV, fPhoton_Corrected_MomY_Col_GeV;

double fQsq_Corrected_GeV, fQsq_Corrected, fW_Corrected, fW_Corrected_GeV, fT_Corrected, fT_Corrected_GeV, fx_Corrected, fy_Corrected, fz_Corrected;

double fWFactor, fA, fFlux_Factor_Col, fFlux_Factor_RF, fJacobian_CM, fJacobian_CM_RF, fJacobian_CM_Col, fZASig_T, fZASig_L, fZASig_LT, fZASig_TT, ftestsig, fZASig_L2;

double fZASigma_UU, fRorySigma_UT, fSigma_Col, fSigma_UUPara, fSig_VR, fSig_L, fSig_T;

double fSig_fpi_6GeV;

double fSigmaPhiS, fSigmaPhi_Minus_PhiS, fSigma2Phi_Minus_PhiS, fSigma3Phi_Minus_PhiS, fSigmaPhi_Plus_PhiS, fSigma2Phi_Plus_PhiS, fSig_Phi_Minus_PhiS, fSig_PhiS, fSig_2Phi_Minus_PhiS, fSig_Phi_Plus_PhiS, fSig_3Phi_Minus_PhiS, fSig_2Phi_Plus_PhiS, fEventWeight, fEventWeightMax, fEventWeightCeil, fEventWeightRn, fZAWFactor, fRR, fPhaseSpaceWeight, fPhaseShiftWeight, fWilliamsWeight, fDedrickWeight, fCatchenWeight, fPhi, fPhiS, fPhi_Corrected, fPhiS_Corrected;

double fElectron_Mom_Sq_RF, fElectron_Mom_Sq_Col, fProton_Mom_Sq_Col, fProton_Mom_Sq_CM, fProton_Mom_Sq_RF, fPhoton_Mom_Sq_Col, fPhoton_Mom_Sq_CM, fPhoton_Mom_Sq_RF, fPion_Mom_Sq_Col, fPion_Mom_Sq_CM, fPion_Mom_Sq_RF, fNeutron_Mom_Sq_Col, fNeutron_Mom_Sq_CM, fNeutron_Mom_Sq_RF, fScatElec_Mom_Sq_Col, fScatElec_Mom_Sq_RF;

double fAsymPhiMinusPhi_S, fAsymPhi_S, fAsym2PhiMinusPhi_S, fAsymPhiPlusPhi_S, fAsym3PhiMinusPhi_S, fAsym2PhiPlusPhi_S;

double fTerm_PhiMinusPhi_S, fTerm_Phi_S, fTerm_2PhiMinusPhi_S, fTerm_PhiPlusPhi_S, fTerm_3PhiMinusPhi_S, fTerm_2PhiPlusPhi_S;

double fAsymPhiMinusPhi_S_Col, fAsymPhi_S_Col, fAsym2PhiMinusPhi_S_Col, fAsymPhiPlusPhi_S_Col, fAsym3PhiMinusPhi_S_Col, fAsym2PhiPlusPhi_S_Col;

double fTerm_PhiMinusPhi_S_Col, fTerm_Phi_S_Col, fTerm_2PhiMinusPhi_S_Col, fTerm_PhiPlusPhi_S_Col, fTerm_3PhiMinusPhi_S_Col, fTerm_2PhiPlusPhi_S_Col;

double fPhi_Pion_LeptonPlane_RF, fCos_Phi_Pion_LeptonPlane_RF, fSin_Phi_Pion_LeptonPlane_RF, fPhi_TargPol_LeptonPlane_RF, fCos_Phi_TargPol_LeptonPlane_RF, fSin_Phi_TargPol_LeptonPlane_RF, fTheta_Pion_Photon_RF, fPhi_Pion_LeptonPlane_Col, fCos_Phi_Pion_LeptonPlane_Col, fSin_Phi_Pion_LeptonPlane_Col, fPhi_TargPol_LeptonPlane_Col, fCos_Phi_TargPol_LeptonPlane_Col, fSin_Phi_TargPol_LeptonPlane_Col, fTheta_Pion_Photon_Col;

double fPhi_Omega_LeptonPlane_RF, fCos_Phi_Omega_LeptonPlane_RF, fSin_Phi_Omega_LeptonPlane_RF, fTheta_Omega_Photon_RF;

double fZASigma_UU_Col, fRorySigma_UT_Col, fSig_Phi_Minus_PhiS_Col, fSig_PhiS_Col, fSig_2Phi_Minus_PhiS_Col, fSig_Phi_Plus_PhiS_Col, fSig_3Phi_Minus_PhiS_Col, fSig_2Phi_Plus_PhiS_Col;

double fepi1, fepi2, fradical;

double fOmega_Energy_CM, fOmega_Mom_CM, fOmega_Energy_CM_GeV, fOmega_Mom_CM_GeV;   

double fMomentum[300];

vector<vector<vector<vector<double>>>> SigPar;

int psf_steps; // Love Preet - Added for phase space factor calculations

double psf_ScatElec_E_Stepsize, psf_ScatElec_Theta_Stepsize, psf_ScatElec_Phi_Stepsize, psf_Ejec_Theta_Stepsize; // Love Preet - Added for phase space factor calculations

double psf_ScatElec_E, psf_ScatElec_Theta, psf_ScatElec_Phi, psf_ScalElec_Mom, psf_Ejectile_Theta, psf_Ejectile_Phi, psf_Q2, psf_W, psf_W2, psf_t, psf_ScatElec_Theta_max, psf_ScatElec_Theta_min, psf_ScatElec_E_max, psf_ScatElec_E_min, psf_Ejectile_Theta_max, psf_Ejectile_Theta_min; // Love Preet - Added for phase space factor calculations

double fScatElec_Energy_Col_max, fScatElec_Energy_Col_min, fScatElec_Theta_Col_max, fScatElec_Theta_Col_min,  f_Ejectile_Theta_Col_max, f_Ejectile_Theta_Col_min, fPSF_org;
// Love Preet - Added for actual phase space factor calculations

double scat_e_px, scat_e_py, scat_e_pz, scat_e_E, ejec_px, ejec_py, ejec_pz, ejec_E, rclH_px, rclH_py, rclH_pz, rclH_E; // Love Preet - Added to be stored in the root tree

double fBeta_Col, fGamma_Col, ftheta_Col; //Love Preet - Added for jacobian calculation in collider frame

double fProb[300] = {    
		     6.03456,    6.02429,    6.01155,    5.99636,    5.97873,    5.95869,    5.93626,    5.91147,    5.88435,    5.85493,
		     5.82325,    5.78935,    5.75326,    5.71504,    5.67472,    5.63235,    5.58799,    5.54169,     5.4935,    5.44347,		   
		     5.39167,    5.33816,    5.28299,    5.22623,    5.16794,    5.10818,    5.04703,    4.98455,    4.92081,    4.85588,
		     4.78982,    4.71692,    4.63621,    4.55583,    4.47582,    4.39621,    4.31702,    4.23828,    4.16002,    4.08227,
		     4.00506,     3.9284,    3.85233,    3.77686,    3.70202,    3.62783,    3.55432,    3.48149,    3.40937,    3.33798,
		     3.26733,    3.19745,    3.12834,    3.06002,    2.99251,    2.92581,    2.85995,    2.79493,    2.73075,    2.66744,
		     2.605,      2.54344,    2.48276,    2.41728,    2.35244,    2.28922,    2.22759,    2.16751,    2.10895,    2.05186,
		     1.99621,    1.94198,    1.88913,    1.83762,    1.78743,    1.73851,    1.69086,    1.64442,    1.59918,    1.55511,			   
		     1.51217,    1.47035,    1.42961,    1.38993,    1.35128,    1.31364,    1.27698,    1.24129,    1.20653,    1.17269,			   
		     1.13973,    1.10765,    1.07642,    1.04601,    1.01626,   0.986934,   0.958443,   0.930759,   0.903861,   0.877727,			   
		     0.852335,   0.827665,   0.803696,   0.780409,   0.757785,   0.735806,   0.714452,   0.693708,   0.673555,   0.653978,			   
		     0.63496,   0.616485,   0.598538,   0.581105,    0.56417,   0.547721,   0.531743,   0.516223,   0.501148,   0.486506,			   
		     0.472284,    0.45847,   0.445054,   0.432024,   0.419368,   0.407077,   0.395125,   0.383513,   0.372237,   0.361289,			   
		     0.350658,   0.340336,   0.330314,   0.320583,   0.311135,   0.301962,   0.293056,   0.284409,   0.276014,   0.267863,			   
		     0.25995,   0.252268,   0.244809,   0.237569,   0.230539,   0.223715,   0.217091,    0.21066,   0.204417,   0.198357,			   
		     0.192474,   0.186763,   0.181219,   0.175838,   0.170615,   0.165545,   0.160623,   0.155843,   0.151185,   0.146666,			   
		     0.142282,   0.138028,   0.133902,   0.129898,   0.126013,   0.122245,   0.118588,   0.115041,   0.111599,   0.108261,			   
		     0.105021,   0.101879,  0.0988298,  0.0958719,  0.0930022,  0.0902182,  0.0875173,   0.084897,  0.0823549,  0.0798886,			   
		     0.0774961,  0.0751749,  0.0729231,  0.0707385,  0.0686192,  0.0665631,  0.0645685,  0.0626335,  0.0607563,  0.0589545,			   
		     0.057219,  0.0555322,   0.053893,  0.0523001,  0.0507522,  0.0492481,  0.0477867,  0.0463667,  0.0449872,  0.0436468,			   
		     0.0423448,  0.0410798,  0.0398511,  0.0386576,  0.0374982,  0.0363722,  0.0352785,  0.0342164,  0.0331849,  0.0321831,			   
		     0.0312104,  0.0302658,  0.0293486,  0.0284581,  0.0275935,   0.026754,  0.0259391,  0.0251479,  0.0243799,  0.0236344,			   
		     0.0229107,  0.0221901,  0.0214923,  0.0208167,  0.0201626,  0.0195293,  0.0189162,  0.0183226,  0.0177478,  0.0171913,			   
		     0.0166525,  0.0161308,  0.0156257,  0.0151366,  0.0146629,  0.0142044,  0.0137603,  0.0133303,  0.0129139,  0.0125107,			   
		     0.0121203,  0.0117421,   0.011376,  0.0110214,   0.010678,  0.0103455,  0.0100234, 0.00971153, 0.00940947, 0.00911693,			   
		     0.00883361,  0.0085592, 0.00829385, 0.00803735, 0.00778883, 0.00754804, 0.00731474,  0.0070887, 0.00686967, 0.00665746,			   
		     0.00645184, 0.00625261, 0.00605957, 0.00587252, 0.00569128, 0.00551567, 0.00534551, 0.00518063, 0.00502086, 0.00486605,			   
		     0.00471604, 0.00457069, 0.00442984, 0.00429336,  0.0041611, 0.00403295, 0.00390877, 0.00378843, 0.00367182, 0.00355882,			   
		     0.00344932, 0.00334321, 0.00324038, 0.00314073, 0.00304505,  0.0029524, 0.00286252, 0.00277533, 0.00269076, 0.00260872,			   
		     0.00252913, 0.00245194, 0.00237706, 0.00230444, 0.00223399, 0.00216566, 0.00209939, 0.00203512, 0.00197277, 0.00191231 };


pim::pim() {
}

pim::pim(int aaa) {

  gen_seed = aaa;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void pim::Initilize() {

  //    fRandom = new TRandom2(0);
  //    fRandom->GetSeed();
  //    fRandom->SetSeed(gen_seed);
	
  fRandom = new TRandom3();

  fRandom->SetSeed(gen_seed);
	
  //	cout << fRandom->GetSeed() << endl;
  //	cout << "Seed Used: " << gen_seed << endl;
	
  //	exit(0);

  allset                                      = false;
  kCalcFermi                                  = false;
  kCalcBremss                                 = false;
  kCalcIon                                    = false;
  kCalcBremssEle                              = false;
  kCalcIonEle                                 = false;
  kFSI                                        = false;
  kMSele                                      = false;
  kMS                                         = false;
  gROOTOut                                    = false; // Set to false by default
  
  // 18/01/23 - The luminosity below is some default assumtpion, more up to date values are set in DEMP prod and depend upon beam energy combinations if they are specified
  // See slide 11 in https://indico.cern.ch/event/1072579/contributions/4796856/attachments/2456676/4210776/CAP-EIC-June-7-2022-Seryi-r2.pdf for more info
  // fLumi                                     = 0.374e33; // Jlab design
  fLumi                                       = 1e33; // 18/01/23, this seems a better default based upon more up to date info, see link above
  fuBcm2                                      = 1.0e-30;
  fPI                                         = 3.1415926;
  fDEG2RAD                                    = fPI/180.0;
  fRAD2DEG                                   = 180.0/fPI;
    
  // SJDK 21/12/22 - Set by .json read in
  //fScatElec_Theta_I                           = 60.0 * fDEG2RAD;
  //fScatElec_Theta_F                           = 175.0 * fDEG2RAD;
  // SJDK 29/11/22 - Updated comment on two variables below
  // Two parameters below are NOT a percentage of the beam energy as previously clamed. This parameter along with _Hi represent the RANGE of enegries over which the scattered electron is generated
  // The range is from 0.5*EBeam to 2.5*EBeam -> Therefore for the phase space calculation, the spread of 2* the incoming beam energy is used in the calculation
  // SJDK 21/12/22 - Set by .json read in
  //fScatElec_E_Lo                              = 0.5;  // NOT a percentage of beam energy
  //fScatElec_E_Hi                              = 2.5;  // NOT a percentage of beam energy
  // Should remove specific pion/omega values here and JUST use EjectileX, should then actually read this in as a parameter, use 0 to 50 as default
  fPion_Theta_I                               = 0.0 * fDEG2RAD;
  fPion_Theta_F                               = 50.0 * fDEG2RAD;
  // SJDK 21/12/22 - Set by .json read in
  //fEjectileX_Theta_I                          = 0.0 * fDEG2RAD;
  //fEjectileX_Theta_F                          = 50.0 * fDEG2RAD;
  fOmega_Theta_I                              = 0.0 * fDEG2RAD; 
  fOmega_Theta_F                              = 360.0 * fDEG2RAD; 
  // 02/06/21 - SJDK
  // Set to 0, now set in PiPlusProd.cc
  fPSF                                     = 0;
  fK                                          = 1000.0;
  fm                                          = 1.0/1000.0;
  fElectron_Mass                              = 0.51099895000;
  fElectron_Mass_GeV                          = fElectron_Mass/1000.0;
  fProton_Mass                                = 938.27208816; // Its is the mass of Proton which in SoLID DVMP is outgoing reocil Proton
  fProton_Mass_GeV                            = fProton_Mass/1000.0;
  fNeutron_Mass                               = 939.5654205; // It is the mass of Neutron. The target is Neutron in SoLID DVMP.
  fNeutron_Mass_GeV                           = fNeutron_Mass/1000.0;
  fRecoilProton_Mass                          = 938.27208816;
  fRecoilProton_Mass_GeV                      = fRecoilProton_Mass/1000.0;
  fPion_Mass                                  = 139.57039 ;
  fPion_Mass_GeV                              = fPion_Mass/1000.0;
  fKaon_Mass                                  = 493.677;
  fKaon_Mass_GeV                              = fKaon_Mass/1000.0;
  fLambda_Mass                                = 1115.683;
  fLambda_Mass_GeV                            = fLambda_Mass/1000.0;
  fSigma_Mass                                 = 1192.642;
  fSigma_Mass_GeV                             = fSigma_Mass/1000.0;
  fOmega_Mass                                 = 782.66;
  fOmega_Mass_GeV                             = fOmega_Mass/1000.0;

  fDiff                                       = 0.00001; // 10/05/23 - Love Preet - Changed from 0.5
    
  // Love Preet - Adding for phase space factor calculations
  psf_steps                                   = 1000.0; //1000.0;
  psf_ScatElec_E_max = std::numeric_limits<double>::min(); // Initialize maxValue for the scattered electron's energy
  psf_ScatElec_E_min = std::numeric_limits<double>::max(); // Initialize minValue for the scattered electron's energy
  psf_ScatElec_Theta_max = std::numeric_limits<double>::min(); // Initialize maxValue for the scattered electron's theta
  psf_ScatElec_Theta_min = std::numeric_limits<double>::max(); // Initialize minValue for the scattered electron's theta
  psf_Ejectile_Theta_max = std::numeric_limits<double>::min(); // Initialize maxValue for the ejectile's theta
  psf_Ejectile_Theta_min = std::numeric_limits<double>::max(); // Initialize minValue for the ejectile's theta
    
  // Love Preet - Added for actual phase space factor calculations
  fScatElec_Energy_Col_max = std::numeric_limits<double>::min(); // Initialize maxValue for the scattered electron's energy
  fScatElec_Energy_Col_min = std::numeric_limits<double>::max(); // Initialize minValue for the scattered electron's energy
  fScatElec_Theta_Col_max = std::numeric_limits<double>::min(); // Initialize maxValue for the scattered electron's theta
  fScatElec_Theta_Col_min = std::numeric_limits<double>::max(); // Initialize minValue for the scattered electron's theta
  f_Ejectile_Theta_Col_max = std::numeric_limits<double>::min(); // Initialize maxValue for the ejectile's theta
  f_Ejectile_Theta_Col_min = std::numeric_limits<double>::max(); // Initialize minValue for the ejectile's theta
  fPSF_org = 0;
    
  // Love Preet - Added to be stored in the root tree
  scat_e_px                                  = 0;
  scat_e_py                                  = 0;
  scat_e_pz                                  = 0;
  scat_e_E                                   = 0;
  ejec_px                                    = 0;
  ejec_py                                    = 0;
  ejec_pz                                    = 0;
  ejec_E                                     = 0;
  rclH_px                                    = 0;
  rclH_py                                    = 0;
  rclH_pz                                    = 0;
  rclH_E                                     = 0;    
  
  //Love Preet - Added for jacobian calculation in collider frame
  fBeta_Col                                  = 0;
  fGamma_Col                                 = 0;
  ftheta_Col                                 = 0;
  // 02/06/21 - SJDK
  // Set to 0, now set in PiPlusProd.cc
  fElectron_Kin_Col_GeV                       = 0;
  fElectron_Kin_Col                           = 0;
  fAlpha                                      = 1./137.036;
  fMom_Ratio                                  = 0.460029;
  fMom_Dif                                    = 0.01;
  fPi                                         = TMath::Pi(); 
  fMandSConserve                              = 0;
  fEnergyConserve                             = 0;
  fXMomConserve                               = 0;
  fYMomConserve                               = 0;
  fZMomConserve                               = 0;
  fXMomConserve_RF                            = 0;
  fYMomConserve_RF                            = 0;
  fZMomConserve_RF                            = 0;
  fEnergyConserve_RF                          = 0;
  fPion_Mom_Same                              = 0;
  fTop_Pion_Mom                               = 0;
  fBot_Pion_Mom                               = 0;
  fS_I_RF                                     = 0;
  fS_F_RF                                     = 0;
  fS_I_Col                                    = 0;
  fS_F_Col                                    = 0;
  fS_I_RF_GeV                                 = 0;
  fS_F_RF_GeV                                 = 0;
  fS_I_Col_GeV                                = 0;
  fS_F_Col_GeV                                = 0;
  fPion_Alpha                                 = 0;
  fPion_Beta                                  = 0;
  fNRecorded                                  = 0;
  fNGenerated                                 = 0;
  fRatio                                      = 0;
  fWLessShell                                 = 0;
  fWLess1P9                                   = 0;
  fWSqNeg                                     = 0;
  fNSigmaNeg                                  = 0;
  fNWeightNeg                                 = 0;
  // SJDK 15/06/21 - Integer counters to check number returning NaN and failing conservation laws added
  fNaN                                        = 0;
  fConserve                                   = 0;
  fNWeightUnphys                              = 0;
  fNWeightReject                              = 0;
  fSolveEvents_0Sol                           = 0;
  fSolveEvents_1Sol                           = 0;
  fSolveEvents_2Sol                           = 0;
  fSDiff                                      = 0;
  fScatElecEnergyLess                         = 0;
  fScatElecThetaLess                          = 0;
  fPionEnergyCMLess                           = 0;
  fSNotEqual                                  = 0;
  fVertex_X                                   = 0;
  fVertex_Y                                   = 0;
  fVertex_Z                                   = 0;
  fProton_Energy_Col                          = 0;
  fProton_Mom_Col                             = 0;
  fProton_Theta_Col                           = 0;
  fProton_Phi_Col                             = 0;
  fProton_MomZ_Col                            = 0;
  fProton_MomX_Col                            = 0;
  fProton_MomY_Col                            = 0;
  fProton_Energy_Col_GeV                      = 0;
  fProton_Mom_Col_GeV                         = 0;
  fProton_MomX_Col_GeV                        = 0;
  fProton_MomY_Col_GeV                        = 0;
  fProton_MomZ_Col_GeV                        = 0;
  fTarget_Energy_Col                          = 0;
  fTarget_Mom_Col                             = 0;
  fTarget_Theta_Col                           = 0;
  fTarget_Phi_Col                             = 0;
  fTarget_MomZ_Col                            = 0;
  fTarget_MomX_Col                            = 0;
  fTarget_MomY_Col                            = 0;
  fTarget_Energy_Col_GeV                      = 0;
  fTarget_Mom_Col_GeV                         = 0;
  fTarget_MomX_Col_GeV                        = 0;
  fTarget_MomY_Col_GeV                        = 0;
  fTarget_MomZ_Col_GeV                        = 0;
  fTarget_Pol0_Col                            = 0;
  fTarget_PolX_Col                            = 0;
  fTarget_PolY_Col                            = 0;
  fTarget_PolZ_Col                            = 0;
  fTarget_Pol0_RF                             = 0;
  fTarget_PolX_RF                             = 0;
  fTarget_PolY_RF                             = 0;
  fTarget_PolZ_RF                             = 0;
  fBetaX_Col_RF                               = 0;
  fBetaY_Col_RF                               = 0;
  fBetaZ_Col_RF                               = 0;
  fBeta_Col_RF                                = 0;
  fGamma_Col_RF                               = 0;
  fProton_MomX_RF                             = 0;
  fProton_MomY_RF                             = 0;
  fProton_MomZ_RF                             = 0;
  fProton_Mom_RF                              = 0;
  fProton_Energy_RF                           = 0;
  fProton_Energy_RF_GeV                       = 0;
  fProton_MomX_RF_GeV                         = 0;
  fProton_MomY_RF_GeV                         = 0;
  fProton_MomZ_RF_GeV                         = 0;
  fProton_Mom_RF_GeV                          = 0;
  fProton_Kin_Col_GeV                         = 0;
  fScatElec_Angle                             = 0;
  fScatElec_Alpha_RF                          = 0;
  fScatElec_Beta_RF                           = 0;
  fRadiation_Lenght_Air                       = 0;
  fElectron_Targ_Thickness                    = 0;
  fElectron_Targ_Thickness_RadLen             = 0;
  fElectron_Targ_BT                           = 0;
  fElectron_Targ_Bremss_Loss                  = 0;
  fElectron_Targ_Ion_Loss                     = 0;
  fElectron_TargWindow_Bremss_Loss            = 0;
  fElectron_TargWindow_Ion_Loss               = 0;
  fElectron_Air_Thickness                     = 0;
  fElectron_Air_Thickness_RadLen              = 0;
  fElectron_Air_BT                            = 0;
  fElectron_Air_Bremss_Loss                   = 0;
  fElectron_Air_Ion_Loss                      = 0;
  fElectron_Corrected_Theta_Col               = 0;
  fElectron_Corrected_Phi_Col                 = 0;
  fElectron_Corrected_Energy_Col              = 0;
  fElectron_Corrected_Mom_Col                 = 0;
  fElectron_Corrected_MomX_Col                = 0;
  fElectron_Corrected_MomY_Col                = 0;
  fElectron_Corrected_MomZ_Col                = 0;
  fElectron_Delta_Mom_Col                     = 0;
  fElectron_Corrected_Energy_Col_GeV          = 0;
  fElectron_Corrected_Mom_Col_GeV             = 0;
  fElectron_Corrected_MomX_Col_GeV            = 0;
  fElectron_Corrected_MomY_Col_GeV            = 0;
  fElectron_Corrected_MomZ_Col_GeV            = 0;
  fElectron_Delta_Mom_Col_GeV                 = 0;
  
  fElectron_Energy_Col                        = 0;
  fElectron_MomZ_Col                          = 0;
  fElectron_MomX_Col                          = 0;
  fElectron_MomY_Col                          = 0;
  fElectron_Theta_Col                         = 0;
  fElectron_Phi_Col                           = 0;
  fElectron_Mom_Col                           = 0;
    
  fElectron_MS_Energy_Col                     = 0;
  fElectron_MS_MomZ_Col                       = 0;
  fElectron_MS_MomX_Col                       = 0;
  fElectron_MS_MomY_Col                       = 0;
  fElectron_MS_Theta_Col                      = 0;
  fElectron_MS_Phi_Col                        = 0;
  fElectron_MS_Mom_Col                        = 0;
  
  fElectron_Energy_Col_GeV                    = 0;
  fElectron_Mom_Col_GeV                       = 0;
  fElectron_MomX_Col_GeV                      = 0;
  fElectron_MomY_Col_GeV                      = 0;
  fElectron_MomZ_Col_GeV                      = 0;
  fScatElec_Targ_Thickness                    = 0;
  fScatElec_Targ_Thickness_RadLen             = 0;
  fScatElec_Targ_BT                           = 0;
  fScatElec_Targ_Bremss_Loss                  = 0;
  fScatElec_Targ_Ion_Loss                     = 0;
  fScatElec_Air_Thickness                     = 0;
  fScatElec_Air_Thickness_RadLen              = 0;
  fScatElec_Air_BT                            = 0;
  fScatElec_Air_Bremss_Loss                   = 0;
  fScatElec_Air_Ion_Loss                      = 0;
  fScatElec_Corrected_Theta_Col               = 0;
  fScatElec_Corrected_Phi_Col                 = 0;
  fScatElec_Corrected_Energy_Col              = 0;
  fScatElec_Corrected_Mom_Col                 = 0;
  fScatElec_Corrected_MomX_Col                = 0;
  fScatElec_Corrected_MomY_Col                = 0;
  fScatElec_Corrected_MomZ_Col                = 0;
  fScatElec_Delta_Mom_Col                     = 0;
  fScatElec_Corrected_Energy_Col_GeV          = 0;
  fScatElec_Corrected_Mom_Col_GeV             = 0;
  fScatElec_Corrected_MomX_Col_GeV            = 0;
  fScatElec_Corrected_MomY_Col_GeV            = 0;
  fScatElec_Corrected_MomZ_Col_GeV            = 0;
  fScatElec_Delta_Mom_Col_GeV                 = 0;
  fScatElec_Energy_Col                        = 0;
  fScatElec_MomZ_Col                          = 0;
  fScatElec_MomX_Col                          = 0;
  fScatElec_MomY_Col                          = 0;
  fScatElec_Theta_Col                         = 0;
  fScatElec_Phi_Col                           = 0;
  fScatElec_Mom_Col                           = 0;
  fScatElec_Energy_Col_GeV                    = 0;
  fScatElec_Mom_Col_GeV                       = 0;
  fScatElec_MomX_Col_GeV                      = 0;
  fScatElec_MomY_Col_GeV                      = 0;
  fScatElec_MomZ_Col_GeV                      = 0;
  
  fScatElec_MS_Energy_Col                     = 0;
  fScatElec_MS_MomZ_Col                       = 0;
  fScatElec_MS_MomX_Col                       = 0;
  fScatElec_MS_MomY_Col                       = 0;
  fScatElec_MS_Theta_Col                      = 0;
  fScatElec_MS_Phi_Col                        = 0;
  fScatElec_MS_Mom_Col                        = 0;
    
  fScatElec_TargWindow_Bremss_Loss            = 0;
  fScatElec_TargWindow_Ion_Loss               = 0;
  fTargWindow_Thickness                       = 0;
  fTargWindow_Thickness_RadLen                = 0;
  fTargWindow_BT                              = 0; 
  fPion_TargWindow_Ion_Loss                   = 0;
  fNeutron_TargWindow_Ion_Loss                = 0;
    
  fPion_MS_Energy_Col                         = 0;
  fPion_MS_MomZ_Col                           = 0;
  fPion_MS_MomX_Col                           = 0;
  fPion_MS_MomY_Col                           = 0;
  fPion_MS_Theta_Col                          = 0;
  fPion_MS_Phi_Col                            = 0;
  fPion_MS_Mom_Col                            = 0;
  
  fPion_Targ_Thickness                        = 0;
  fPion_Targ_Thickness_RadLen                 = 0;
  fPion_Targ_BT                               = 0;
  fPion_Targ_Bremss_Loss                      = 0;
  fPion_Targ_Ion_Loss                         = 0;
  fPion_Air_Thickness                         = 0;
  fPion_Air_Thickness_RadLen                  = 0;
  fPion_Air_BT                                = 0;
  fPion_Air_Bremss_Loss                       = 0;
  fPion_Air_Ion_Loss                          = 0;
  fPion_Corrected_Theta_Col                   = 0;
  fPion_Corrected_Phi_Col                     = 0;
  fPion_Corrected_Energy_Col                  = 0;
  fPion_Corrected_Mom_Col                     = 0;
  fPion_Corrected_MomX_Col                    = 0;
  fPion_Corrected_MomY_Col                    = 0;
  fPion_Corrected_MomZ_Col                    = 0;
  fPion_Delta_Mom_Col                         = 0;
  fPion_Corrected_Energy_Col_GeV              = 0;
  fPion_Corrected_Mom_Col_GeV                 = 0;
  fPion_Corrected_MomX_Col_GeV                = 0;
  fPion_Corrected_MomY_Col_GeV                = 0;
  fPion_Corrected_MomZ_Col_GeV                = 0;
  fPion_Delta_Mom_Col_GeV                     = 0;
    
  fPion_Energy_Col                            = 0;
  fPion_MomZ_Col                              = 0;
  fPion_MomX_Col                              = 0;
  fPion_MomY_Col                              = 0;
  fPion_Theta_Col                             = 0;
  fPion_Phi_Col                               = 0;
  fPion_Mom_Col                               = 0;
  fPion_Energy_Col_GeV                        = 0;
  fPion_Mom_Col_GeV                           = 0;
  fPion_MomX_Col_GeV                          = 0;
  fPion_MomY_Col_GeV                          = 0;
  fPion_MomZ_Col_GeV                          = 0;
    
  fPion_FSI_Energy_Col                        = 0;
  fPion_FSI_MomZ_Col                          = 0;
  fPion_FSI_MomX_Col                          = 0;
  fPion_FSI_MomY_Col                          = 0;
  fPion_FSI_Theta_Col                         = 0;
  fPion_FSI_Phi_Col                           = 0;
  fPion_FSI_Mom_Col                           = 0;
  fPion_FSI_Energy_Col_GeV                    = 0;
  fPion_FSI_Mom_Col_GeV                       = 0;
  fPion_FSI_MomX_Col_GeV                      = 0;
  fPion_FSI_MomY_Col_GeV                      = 0;
  fPion_FSI_MomZ_Col_GeV                      = 0;
  
  fNeutron_MS_Energy_Col                      = 0;
  fNeutron_MS_MomZ_Col                        = 0;
  fNeutron_MS_MomX_Col                        = 0;
  fNeutron_MS_MomY_Col                        = 0;
  fNeutron_MS_Theta_Col                       = 0;
  fNeutron_MS_Phi_Col                         = 0;
  fNeutron_MS_Mom_Col                         = 0;

  fKaon_Energy_Col                            = 0;       
  fKaon_MomZ_Col                              = 0;       
  fKaon_MomX_Col                              = 0;       
  fKaon_MomY_Col                              = 0;       
  fKaon_Theta_Col                             = 0;       
  fKaon_Phi_Col                               = 0;       
  fKaon_Mom_Col                               = 0;       
  fKaon_Energy_Col_GeV                        = 0;       
  fKaon_Mom_Col_GeV                           = 0;       
  fKaon_MomX_Col_GeV                          = 0;       
  fKaon_MomY_Col_GeV                          = 0;       
  fKaon_MomZ_Col_GeV                          = 0;
  fScathad_Energy_Col                         = 0;    
  fScathad_MomZ_Col                           = 0;    
  fScathad_MomX_Col                           = 0;    
  fScathad_MomY_Col                           = 0;    
  fScathad_Theta_Col                          = 0;    
  fScathad_Phi_Col                            = 0;    
  fScathad_Mom_Col                            = 0;    
  fScathad_Energy_Col_GeV                     = 0;    
  fScathad_Mom_Col_GeV                        = 0;    
  fScathad_MomX_Col_GeV                       = 0;    
  fScathad_MomY_Col_GeV                       = 0;    
  fScathad_MomZ_Col_GeV                       = 0;
  
  fNeutron_Targ_Thickness                     = 0;
  fNeutron_Targ_Thickness_RadLen              = 0;
  fNeutron_Targ_BT                            = 0;
  fNeutron_Targ_Bremss_Loss                   = 0;
  fNeutron_Targ_Ion_Loss                      = 0;
  fNeutron_Air_Thickness                      = 0;
  fNeutron_Air_Thickness_RadLen               = 0;
  fNeutron_Air_BT                             = 0;
  fNeutron_Air_Bremss_Loss                    = 0;
  fNeutron_Air_Ion_Loss                       = 0;
  fNeutron_Corrected_Theta_Col                = 0;
  fNeutron_Corrected_Phi_Col                  = 0;
  fNeutron_Corrected_Energy_Col               = 0;
  fNeutron_Corrected_Mom_Col                  = 0;
  fNeutron_Corrected_MomX_Col                 = 0;
  fNeutron_Corrected_MomY_Col                 = 0;
  fNeutron_Corrected_MomZ_Col                 = 0;
  fNeutron_Delta_Mom_Col                      = 0;
  fNeutron_Corrected_Energy_Col_GeV           = 0;
  fNeutron_Corrected_Mom_Col_GeV              = 0;
  fNeutron_Corrected_MomX_Col_GeV             = 0;
  fNeutron_Corrected_MomY_Col_GeV             = 0;
  fNeutron_Corrected_MomZ_Col_GeV             = 0;
  fNeutron_Delta_Mom_Col_GeV                  = 0;
  fNeutron_Energy_Col                         = 0;
  fNeutron_MomZ_Col                           = 0;
  fNeutron_MomX_Col                           = 0;
  fNeutron_MomY_Col                           = 0;
  fNeutron_Theta_Col                          = 0;
  fNeutron_Phi_Col                            = 0;
  fNeutron_Mom_Col                            = 0;
  fNeutron_Energy_Col_GeV                     = 0;
  fNeutron_Mom_Col_GeV                        = 0;
  fNeutron_MomX_Col_GeV                       = 0;
  fNeutron_MomY_Col_GeV                       = 0;
  fNeutron_MomZ_Col_GeV                       = 0;
  fRecoilProton_Targ_Thickness                = 0;
  fRecoilProton_Targ_Thickness_RadLen         = 0;
  fRecoilProton_Targ_BT                       = 0;
  fRecoilProton_Targ_Bremss_Loss              = 0;
  fRecoilProton_Targ_Ion_Loss                 = 0;
  fRecoilProton_Air_Thickness                 = 0;
  fRecoilProton_Air_Thickness_RadLen          = 0;
  fRecoilProton_Air_BT                        = 0;
  fRecoilProton_Air_Bremss_Loss               = 0;
  fRecoilProton_Air_Ion_Loss                  = 0;
  fRecoilProton_Theta_Col                     = 0;
  fRecoilProton_Phi_Col                       = 0;
  fRecoilProton_Energy_Col                    = 0;
  fRecoilProton_Mom_Col                       = 0;
  fRecoilProton_MomX_Col                      = 0;
  fRecoilProton_MomY_Col                      = 0;
  fRecoilProton_MomZ_Col                      = 0;
  fRecoilProton_Corrected_Energy_Col          = 0;
  fRecoilProton_Corrected_Mom_Col             = 0;
  fRecoilProton_Corrected_MomX_Col            = 0;
  fRecoilProton_Corrected_MomY_Col            = 0;
  fRecoilProton_Corrected_MomZ_Col            = 0;
  fRecoilProton_Corrected_Theta_Col           = 0;
  fRecoilProton_Corrected_Phi_Col             = 0;
  fRecoilProton_Delta_Mom_Col                 = 0;
  fRecoilProton_Energy_Col_GeV                = 0;
  fRecoilProton_Mom_Col_GeV                   = 0;
  fRecoilProton_MomX_Col_GeV                  = 0;
  fRecoilProton_MomY_Col_GeV                  = 0;
  fRecoilProton_MomZ_Col_GeV                  = 0;
  fRecoilProton_Corrected_Energy_Col_GeV      = 0;
  fRecoilProton_Corrected_Mom_Col_GeV         = 0;
  fRecoilProton_Corrected_MomX_Col_GeV        = 0;
  fRecoilProton_Corrected_MomY_Col_GeV        = 0;
  fRecoilProton_Corrected_MomZ_Col_GeV        = 0;
  fRecoilProton_Delta_Mom_Col_GeV             = 0;
  fSSAsym                                     = 0;
  fSineAsym                                   = 0;
  fT                                          = 0;
  fT_GeV                                      = 0;
  fProton_Kin_Col                             = 0;
  fQsq_Value                                  = 0;
  fQsq_Dif                                    = 0;
  fQsq_GeV                                    = 0;
  fQsq                                        = 0;
  fW_GeV_Col                                  = 0;
  fW_Col                                      = 0;
  fW                                          = 0;   
  fW_GeV                                      = 0;
  fW                                          = 0;   
  fW_GeV                                      = 0;
  fW_Prime_GeV                                = 0;
  fW_Corrected_Prime_GeV                      = 0;
  fWSq                                        = 0;   
  fWSq_GeV                                    = 0;
  fWSq_PiN                                    = 0;   
  fWSq_PiN_GeV                                = 0;
  fWSq_Top_PiN_GeV                            = 0;
  fWSq_Bot_PiN_GeV                            = 0;
  fScatElec_Mom_RF                            = 0;
  fScatElec_Mom_RF_GeV                        = 0;
  fScatElec_Energy_RF                         = 0;
  fScatElec_Energy_RF_GeV                     = 0;
  fElec_ScatElec_Theta_RF                     = 0;
  fScatElec_Cone_Phi_RF                       = 0;
  fScatElec_Theta_RF                          = 0;
  fScatElec_Phi_RF                            = 0;
  fScatElec_MomX_RF                           = 0;
  fScatElec_MomZ_RF                           = 0;
  fScatElec_MomY_RF                           = 0;
  fElectron_Energy_RF                         = 0;
  fElectron_Mom_RF                            = 0;
  fElectron_Theta_RF                          = 0;
  fElectron_Phi_RF                            = 0;
  fElectron_MomX_RF                           = 0;
  fElectron_MomZ_RF                           = 0;
  fElectron_MomY_RF                           = 0;
  fPhoton_Energy_RF                           = 0;
  fPhoton_Mom_RF                              = 0;
  fPhoton_Energy_RF_GeV                       = 0;
  fPhoton_Mom_RF_GeV                          = 0;
  fProton_Energy_CM                           = 0;
  fProton_Energy_CM_GeV                       = 0;
  fProton_Mom_CM                              = 0;
  fProton_Mom_CM_GeV                          = 0;
  fPhoton_Energy_CM                           = 0;
  fPhoton_Mom_CM                              = 0;
  fPhoton_Energy_CM_GeV                       = 0;
  fPhoton_Mom_CM_GeV                          = 0;
  fPion_Theta_CM                              = 0;
  fPion_Phi_CM                                = 0;
  fPion_Energy_CM                             = 0;
  fPion_Mom_CM                                = 0;
  fPion_Energy_CM_GeV                         = 0;
  fPion_Mom_CM_GeV                            = 0;
  fBeta_CM_RF                                 = 0;
  fGamma_CM_RF                                = 0;
  fNeutron_Theta_CM                           = 0;
  fNeutron_Phi_CM                             = 0;
  fNeutron_Energy_CM                          = 0;
  fNeutron_Energy_CM_GeV                      = 0;
  fNeutron_Mom_CM                             = 0;
  fNeutron_Mom_CM_GeV                         = 0;
  fPhoton_MomZ_RF                             = 0;
  fPhoton_MomX_RF                             = 0;
  fPhoton_MomY_RF                             = 0;
  fPhoton_Theta_RF                            = 0;
  fPhoton_Phi_RF                              = 0;
  fPion_Energy_RF                             = 0;
  fPion_Mom_RF                                = 0;
  fPion_Energy_RF_GeV                         = 0;
  fPion_Mom_RF_GeV                            = 0;
  fPiqVec_Theta_RF                            = 0;
  fPion_Mom_RF                                = 0;
  fPion_Mom_RF_GeV                            = 0;
  fPion_MomX_RF                               = 0;
  fPion_MomY_RF                               = 0;
  fPion_MomZ_RF                               = 0;
  fPion_Theta_RF                              = 0;
  fPion_Phi_RF                                = 0;
  fPion_MomX_RF_GeV                           = 0;
  fPion_MomY_RF_GeV                           = 0;
  fPion_MomZ_RF_GeV                           = 0;
  fT_Para                                     = 0;
  fT_Para_GeV                                 = 0;
  fEpsilon                                    = 0;
  fx                                          = 0;
  fy                                          = 0;
  fz                                          = 0;
  fNeutron_Energy_RF                          = 0;
  fNeutron_Energy_RF_GeV                      = 0;
  fNeutron_Mom_RF                             = 0;
  fNeutron_Mom_RF_GeV                         = 0;
  fNeutron_qVec_Theta_RF                      = 0;
  fNeutron_MomX_RF                            = 0;
  fNeutron_MomY_RF                            = 0;
  fNeutron_MomZ_RF                            = 0;
  fNeutron_Theta_RF                           = 0;
  fNeutron_Phi_RF                             = 0;
  fRecoilProton_Energy_RF                     = 0;
  fRecoilProton_Mom_RF                        = 0;
  fRecoilProton_MomX_RF                       = 0;
  fRecoilProton_MomY_RF                       = 0;
  fRecoilProton_MomZ_RF                       = 0;
  fRecoilProton_Energy_RF_GeV                 = 0;
  fRecoilProton_Mom_RF_GeV                    = 0;
  fRecoilProton_MomX_RF_GeV                   = 0;
  fRecoilProton_MomY_RF_GeV                   = 0;
  fRecoilProton_MomZ_RF_GeV                   = 0;
  fRecoilProton_Theta_RF                      = 0;
  fRecoilProton_Phi_RF                        = 0;
  fElectron_MomX_RF_GeV                       = 0;
  fElectron_MomY_RF_GeV                       = 0;
  fElectron_MomZ_RF_GeV                       = 0;
  fPhoton_MomX_RF_GeV                         = 0;
  fPhoton_MomY_RF_GeV                         = 0;
  fPhoton_MomZ_RF_GeV                         = 0;
  fScatElec_MomX_RF_GeV                       = 0;
  fScatElec_MomY_RF_GeV                       = 0;
  fScatElec_MomZ_RF_GeV                       = 0;
  fNeutron_MomX_RF_GeV                        = 0;
  fNeutron_MomY_RF_GeV                        = 0;
  fNeutron_MomZ_RF_GeV                        = 0;
  fPhoton_MomX_Col_GeV                        = 0;
  fPhoton_MomY_Col_GeV                        = 0;
  fPhoton_MomZ_Col_GeV                        = 0;
  fPion_MomX_Col_GeV                          = 0;
  fPion_MomY_Col_GeV                          = 0;
  fPion_MomZ_Col_GeV                          = 0; 
  fPhoton_Theta_Col                           = 0;
  fPhoton_Phi_Col                             = 0;
  fPhoton_Energy_Col                          = 0;
  fPhoton_Mom_Col                             = 0;
  fPhoton_MomX_Col                            = 0;
  fPhoton_MomZ_Col                            = 0;
  fPhoton_MomY_Col                            = 0;
  fPhoton_Energy_Col_GeV                      = 0;
  fPhoton_Mom_Col_GeV                         = 0;
  fPhoton_MomX_Col_GeV                        = 0;
  fPhoton_MomZ_Col_GeV                        = 0;
  fPhoton_MomY_Col_GeV                        = 0;
  fWFactor                                    = 0;
  fA                                          = 0;
  fZASigma_UU                                 = 0;
  fRorySigma_UT                               = 0;
  fSigma_Col                                  = 0;
  fSigma_UUPara                               = 0;
  fSig_VR                                     = 0;

  fSig_fpi_6GeV                               = 0;

  fSig_L                                      = 0;
  fSig_T                                      = 0;
  fSigmaPhiS                                  = 0;
  fSigmaPhi_Minus_PhiS                        = 0;
  fSigma2Phi_Minus_PhiS                       = 0;
  fSigma3Phi_Minus_PhiS                       = 0;
  fSigmaPhi_Plus_PhiS                         = 0;
  fSigma2Phi_Plus_PhiS                        = 0;
  fSig_Phi_Minus_PhiS                         = 0;
  fSig_PhiS                                   = 0;
  fSig_2Phi_Minus_PhiS                        = 0;
  fSig_Phi_Plus_PhiS                          = 0;
  fSig_3Phi_Minus_PhiS                        = 0;
  fSig_2Phi_Plus_PhiS                         = 0;
  fEventWeight                                = 0;
  fEventWeightMax                             = 0;
  fEventWeightCeil                            = 0; // SJDK 11/05/21 - This is the maximum value found with the old method that is used to get the new unit weight
  fEventWeightRn                              = 0; // SJDK 11/05/21 -Random number to compare determined weight to
  fWilliamsWeight                             = 0;
  fDedrickWeight                              = 0;
  fCatchenWeight                              = 0;
  fFlux_Factor_Col                            = 0;
  fFlux_Factor_RF                             = 0;
  fJacobian_CM                                = 0;
  fJacobian_CM_RF                             = 0;
  fJacobian_CM_Col                            = 0;
  fZASig_T                                    = 0;
  fZASig_L                                    = 0;
  fZASig_L2                                   = 0;
  fZASig_LT                                   = 0;
  fZASig_TT                                   = 0;
  fPhi                                        = 0;
  fPhiS                                       = 0;
  fPhi_Corrected                              = 0;
  fPhiS_Corrected                             = 0;      
  
  fQsq_Corrected_GeV                          = 0;      
  fQsq_Corrected                              = 0;      
  fW_Corrected                                = 0;      
  fW_Corrected_GeV                            = 0;      
  fT_Corrected                                = 0;      
  fT_Corrected_GeV                            = 0;      
  fx_Corrected                                = 0;      
  fy_Corrected                                = 0;      
  fz_Corrected                                = 0;      
  
  fAsymPhiMinusPhi_S                          = 0;
  fAsymPhi_S                                  = 0;
  fAsym2PhiMinusPhi_S                         = 0;
  fAsymPhiPlusPhi_S                           = 0;
  fAsym3PhiMinusPhi_S                         = 0;
  fAsym2PhiPlusPhi_S                          = 0;
    
  fAsymPhiMinusPhi_S_Col                      = 0;
  fAsymPhi_S_Col                              = 0;
  fAsym2PhiMinusPhi_S_Col                     = 0;
  fAsymPhiPlusPhi_S_Col                       = 0;
  fAsym3PhiMinusPhi_S_Col                     = 0;
  fAsym2PhiPlusPhi_S_Col                      = 0;
  
  fTerm_PhiMinusPhi_S                         = 0;
  fTerm_Phi_S                                 = 0;
  fTerm_2PhiMinusPhi_S                        = 0;
  fTerm_PhiPlusPhi_S                          = 0;
  fTerm_3PhiMinusPhi_S                        = 0;
  fTerm_2PhiPlusPhi_S                         = 0;
    
  fTerm_PhiMinusPhi_S_Col                     = 0;
  fTerm_Phi_S_Col                             = 0;
  fTerm_2PhiMinusPhi_S_Col                    = 0;
  fTerm_PhiPlusPhi_S_Col                      = 0;
  fTerm_3PhiMinusPhi_S_Col                    = 0;
  fTerm_2PhiPlusPhi_S_Col                     = 0;
  
  fPhoton_Corrected_Theta_Col                 = 0;
  fPhoton_Corrected_Phi_Col                   = 0;
  fPhoton_Corrected_Energy_Col                = 0;
  fPhoton_Corrected_Mom_Col                   = 0;
  fPhoton_Corrected_MomX_Col                  = 0;
  fPhoton_Corrected_MomZ_Col                  = 0;
  fPhoton_Corrected_MomY_Col                  = 0;
  fPhoton_Corrected_Energy_Col_GeV            = 0;
  fPhoton_Corrected_Mom_Col_GeV               = 0;
  fPhoton_Corrected_MomX_Col_GeV              = 0;
  fPhoton_Corrected_MomZ_Col_GeV              = 0;
  fPhoton_Corrected_MomY_Col_GeV              = 0;
  fPhi_Pion_LeptonPlane_RF                    = 0;
  fCos_Phi_Pion_LeptonPlane_RF                = 0;
  fSin_Phi_Pion_LeptonPlane_RF                = 0;

  fPhi_Omega_LeptonPlane_RF                   = 0;
  fCos_Phi_Omega_LeptonPlane_RF               = 0;
  fSin_Phi_Omega_LeptonPlane_RF               = 0;
  fTheta_Omega_Photon_RF                      = 0;

  fOmega_Energy_CM                            = 0;
  fOmega_Mom_CM                               = 0;
  fOmega_Energy_CM_GeV                        = 0;
  fOmega_Mom_CM_GeV                           = 0;

  fPhi_TargPol_LeptonPlane_RF                 = 0;
  fCos_Phi_TargPol_LeptonPlane_RF             = 0;
  fSin_Phi_TargPol_LeptonPlane_RF             = 0;
  fTheta_Pion_Photon_RF                       = 0;
  fPhi_Pion_LeptonPlane_Col                   = 0;
  fCos_Phi_Pion_LeptonPlane_Col               = 0;
  fSin_Phi_Pion_LeptonPlane_Col               = 0;
  fPhi_TargPol_LeptonPlane_Col                = 0;
  fCos_Phi_TargPol_LeptonPlane_Col            = 0;
  fSin_Phi_TargPol_LeptonPlane_Col            = 0;
  fTheta_Pion_Photon_Col                      = 0;
  fZASigma_UU_Col                             = 0;
  fRorySigma_UT_Col                           = 0;
  fSig_Phi_Minus_PhiS_Col                     = 0;
  fSig_PhiS_Col                               = 0;
  fSig_2Phi_Minus_PhiS_Col                    = 0;
  fSig_Phi_Plus_PhiS_Col                      = 0;
  fSig_3Phi_Minus_PhiS_Col                    = 0;
  fSig_2Phi_Plus_PhiS_Col                     = 0;
  // SJDK 08/02/22 - New variables Ali added for conservation law checks
  conserve                                    = 0; // This is the number that PASS both conservation check
  ene                                         = 0; // The number that FAIL due to failing energy conservation check ONLY
  mom                                         = 0; // The number that FAIL due to failing the momentum conservation check ONLY
  ene_mom                                     = 0; // The number that FAIL BOTH energy and momentum check
  mom_px                                      = 0; // Fail due to px only
  mom_py                                      = 0; // Fail due to py only
  mom_pz                                      = 0; // Fail due to pz only
  mom_pxpy                                    = 0; // Fail due to px and py
  mom_pxpz                                    = 0; // Fail due to px and pz
  mom_pypz                                    = 0; // Fail due to py and pz
  mom_pxpypz                                  = 0; // Fail due to px, py and pz

}

//---------------------------------------------------------
double pim::fermiMomentum() {

  double fMom;
  bool kFermi = true;
  while ( kFermi ) {
    double fProton_Rand_Mom_Col      = fRandom->Uniform( 0, 300.0);
    double fProton_Rand_Mom_Col_Prob = fRandom->Uniform( fProb[299], fProb[0] );
    int    fProton_Mom_Int           = std::ceil( fProton_Rand_Mom_Col );
    double f3He_Value                = fProb[ fProton_Mom_Int - 1 ];
  
    if ( fProton_Rand_Mom_Col_Prob <= f3He_Value ) {
      fMom = fProton_Rand_Mom_Col;	
      kFermi = false;
    }
  }

  cout << "Fermi momentum check: " << fMom << endl; 
  
  return fMom;
}

///*--------------------------------------------------*/ 
//-> 10/05/23 - Love added a slimmed down, simpler to read version of the CheckLaws fn
// 
// To check the conservation of the energy and momentum, there two methods avalaible:
// Method 1: Give the four-vectors of the initial and final states partciles, 
//           tolerance factor will be defaulted 1e-6 MeV
// Method 2: Give the four-vectors of the initial and final states partciles, 
//           and the prefered tolerance factor.
//
int pim::CheckLaws(TLorentzVector P_E0, TLorentzVector P_t, TLorentzVector P_e, TLorentzVector P_pim, TLorentzVector P_pro) {

  double energy_check = (P_t.E() + P_E0.E()) - (P_e.E()+P_pim.E()+P_pro.E());
  double px_check =(P_t.Px() + P_E0.Px()) - (P_e.Px()+P_pim.Px()+P_pro.Px());
  double py_check =(P_t.Py() + P_E0.Py()) - (P_e.Py()+P_pim.Py()+P_pro.Py());
  double pz_check =(P_t.Pz() + P_E0.Pz()) - (P_e.Pz()+P_pim.Pz()+P_pro.Pz());
  
  Int_t err = -1;
  if( fabs( energy_check ) < fDiff  && fabs( px_check ) < fDiff  &&  fabs( py_check ) < fDiff && fabs( pz_check ) < fDiff){  // if both momentum components and energy pass the conservation check (simultaneously)
    conserve++;
    err = 1;
  }
  
  else{
    if( (fabs( px_check ) > fDiff) || (fabs( py_check ) > fDiff) || (fabs( pz_check ) > fDiff) ){ // If momentum check fails, check if energy also failed, add counters accordingly  	

      // Check components
      if( fabs( px_check ) > fDiff){
	// px failed, check py
	if( fabs( py_check ) > fDiff){
	  // py failed, check pz
	  if( fabs( pz_check ) > fDiff){
	    // pz failed, all 3 failed
	    mom_pxpypz++; // All failed
	  }
	  else mom_pxpy++; // px and py failed
	}
	else if( fabs( py_check ) < fDiff){
	  // py passed, check pz
	  if( fabs( pz_check ) > fDiff){
	    mom_pxpz++; // px and pz failed
	  }
	  else mom_px++; // px failed, py and pz passd
	}
      }
     
      else if ( abs( px_check ) < fDiff){
	// px passed, check py
	if( fabs( py_check ) > fDiff){
	  // py failed, check pz
	  if( fabs( pz_check ) > fDiff){
	    mom_pypz++; // py and pz failed
	  }
	  else mom_py++; // Just py failed
	}
	else mom_pz++; // Only option left (since we know one of them failed) is that pz failed
      }
      
      if( fabs( energy_check ) > fDiff){
	ene_mom++; // Both failed
      }
      else if( abs( energy_check ) < fDiff){
	mom++; // Only momentum failed
      }
    }    
    else ene++; // If check did not pass, but it wasn't momentum, must have been energy, add to counter
  }
  
  return err;
}

int pim::CheckLaws(TLorentzVector P_E0, TLorentzVector P_t, TLorentzVector P_e, TLorentzVector P_pim, TLorentzVector P_pro, double fDiff_E) {

  double energy_check = (P_t.E() + P_E0.E()) - (P_e.E()+P_pim.E()+P_pro.E());
  double px_check =(P_t.Px() + P_E0.Px()) - (P_e.Px()+P_pim.Px()+P_pro.Px());
  double py_check =(P_t.Py() + P_E0.Py()) - (P_e.Py()+P_pim.Py()+P_pro.Py());
  double pz_check =(P_t.Pz() + P_E0.Pz()) - (P_e.Pz()+P_pim.Pz()+P_pro.Pz());
                                                                                                 
  Int_t err = -1;
  if( fabs( energy_check ) < fDiff_E  && fabs( px_check ) < fDiff_E  &&  fabs( py_check ) < fDiff_E && fabs( pz_check ) < fDiff_E){  // if both momentum components and energy pass the conservation check (simultaneously)
    conserve++;
    err = 1;
  }

  else{
    if( (fabs( px_check ) > fDiff_E) || (fabs( py_check ) > fDiff_E) || (fabs( pz_check ) > fDiff_E) ){ // If momentum check fails, check if energy also failed, add counters accordingly  	

      // Check components
      if( fabs( px_check ) > fDiff_E){
	// px failed, check py
	if( fabs( py_check ) > fDiff_E){
	  // py failed, check pz
	  if( fabs( pz_check ) > fDiff_E){
	    // pz failed, all 3 failed
	    mom_pxpypz++; // All failed
	  }
	  else mom_pxpy++; // px and py failed
	}
	else if( fabs( py_check ) < fDiff_E){
	  // py passed, check pz
	  if( fabs( pz_check ) > fDiff_E){
	    mom_pxpz++; // px and pz failed
	  }
	  else mom_px++; // px failed, py and pz passd
	}
      }
     
      else if ( abs( px_check ) < fDiff_E){
	// px passed, check py
	if( fabs( py_check ) > fDiff_E){
	  // py failed, check pz
	  if( fabs( pz_check ) > fDiff_E){
	    mom_pypz++; // py and pz failed
	  }
	  else mom_py++; // Just py failed
	}
	else mom_pz++; // Only option left (since we know one of them failed) is that pz failed
      }
      
      if( fabs( energy_check ) > fDiff_E){
	ene_mom++; // Both failed
      }
      else if( abs( energy_check ) < fDiff_E){
	mom++; // Only momentum failed
      }
    }    
    else ene++; // If check did not pass, but it wasn't momentum, must have been energy, add to counter
  }
  
  return err;
}

///****************************************

void pim::setrootfile( string rootFile ){ 


  ///****************************************
  /// Bill: re-delcreation of f1 is fixed,
  ///       all object function calls are switched pointer function calls

  f = new TFile(rootFile.c_str(),"recreate"); 

  t1 = new TTree();
  t1->SetName("t1");

  // 6 Particles, electron, proton, scat. electron, photon, pion and neutron in col frame
  //-----------------------------------------------------------------------------------------------------

  /* t1->Branch("Proton_Theta_Col",                          &fProton_Theta_Col,                          "fProton_Theta_Col/D"); */
  /* t1->Branch("Proton_Phi_Col",                            &fProton_Phi_Col,                            "fProton_Phi_Col/D"); */
  /* t1->Branch("Proton_Energy_Col_GeV",                     &fProton_Energy_Col_GeV,                     "fProton_Energy_Col_GeV/D"); */
  /* t1->Branch("Proton_Mom_Col_GeV",                        &fProton_Mom_Col_GeV,                        "fProton_Mom_Col_GeV/D"); */
  /* t1->Branch("Proton_MomZ_Col_GeV",                       &fProton_MomZ_Col_GeV,                       "fProton_MomZ_Col_GeV/D"); */
  /* t1->Branch("Proton_MomX_Col_GeV",                       &fProton_MomX_Col_GeV,                       "fProton_MomX_Col_GeV/D"); */
  /* t1->Branch("Proton_MomY_Col_GeV",                       &fProton_MomY_Col_GeV,                       "fProton_MomY_Col_GeV/D"); */

  /* t1->Branch("Electron_Theta_Col",                        &fElectron_Theta_Col,                        "fElectron_Theta_Col/D"); */
  /* t1->Branch("Electron_Phi_Col",                          &fElectron_Phi_Col,                          "fElectron_Phi_Col/D"); */
  /* t1->Branch("Electron_Energy_Col_GeV",                   &fElectron_Energy_Col_GeV,                   "fElectron_Energy_Col_GeV/D"); */
  /* t1->Branch("Electron_Mom_Col_GeV",                      &fElectron_Mom_Col_GeV,                      "fElectron_Mom_Col_GeV/D"); */
  /* t1->Branch("Electron_MomX_Col_GeV",                     &fElectron_MomX_Col_GeV,                     "fElectron_MomX_Col_GeV/D"); */
  /* t1->Branch("Electron_MomY_Col_GeV",                     &fElectron_MomY_Col_GeV,                     "fElectron_MomY_Col_GeV/D"); */
  /* t1->Branch("Electron_MomZ_Col_GeV",                     &fElectron_MomZ_Col_GeV,                     "fElectron_MomZ_Col_GeV/D"); */

  t1->Branch("ScatElec_Theta_Col",                        &fScatElec_Theta_Col,                        "fScatElec_Theta_Col/D");
  t1->Branch("ScatElec_Phi_Col",                          &fScatElec_Phi_Col,                          "fScatElec_Phi_Col/D");
  t1->Branch("ScatElec_Energy_Col_GeV",                   &fScatElec_Energy_Col_GeV,                   "fScatElec_Energy_Col_GeV/D");
  t1->Branch("ScatElec_Mom_Col_GeV",                      &fScatElec_Mom_Col_GeV,                      "fScatElec_Mom_Col_GeV/D");
  t1->Branch("ScatElec_MomX_Col_GeV",                     &fScatElec_MomX_Col_GeV,                     "fScatElec_MomX_Col_GeV/D");
  t1->Branch("ScatElec_MomY_Col_GeV",                     &fScatElec_MomY_Col_GeV,                     "fScatElec_MomY_Col_GeV/D");
  t1->Branch("ScatElec_MomZ_Col_GeV",                     &fScatElec_MomZ_Col_GeV,                     "fScatElec_MomZ_Col_GeV/D");

  /* t1->Branch("Photon_Theta_Col",                          &fPhoton_Theta_Col,                          "fPhoton_Theta_Col/D"); */
  /* t1->Branch("Photon_Phi_Col",                            &fPhoton_Phi_Col,                            "fPhoton_Phi_Col/D"); */
  /* t1->Branch("Photon_Energy_Col_GeV",                     &fPhoton_Energy_Col_GeV,                     "fPhoton_Energy_Col_GeV/D"); */
  /* t1->Branch("Photon_Mom_Col_GeV",                        &fPhoton_Mom_Col_GeV,                        "fPhoton_Mom_Col_GeV/D"); */
  /* t1->Branch("Photon_MomX_Col_GeV",                       &fPhoton_MomX_Col_GeV,                       "fPhoton_MomX_Col_GeV/D"); */
  /* t1->Branch("Photon_MomY_Col_GeV",                       &fPhoton_MomY_Col_GeV,                       "fPhoton_MomY_Col_GeV/D"); */
  /* t1->Branch("Photon_MomZ_Col_GeV",                       &fPhoton_MomZ_Col_GeV,                       "fPhoton_MomZ_Col_GeV/D"); */

  t1->Branch("Pion_Theta_Col",                            &fPion_Theta_Col,                            "fPion_Theta_Col/D");
  t1->Branch("Pion_Phi_Col",                              &fPion_Phi_Col,                              "fPion_Phi_Col/D");
  t1->Branch("Pion_Energy_Col_GeV",                       &fPion_Energy_Col_GeV,                       "fPion_Energy_Col_GeV/D");
  t1->Branch("Pion_Mom_Col_GeV",                          &fPion_Mom_Col_GeV,                          "fPion_Mom_Col_GeV/D");
  t1->Branch("Pion_MomX_Col_GeV",                         &fPion_MomX_Col_GeV,                         "fPion_MomX_Col_GeV/D");
  t1->Branch("Pion_MomY_Col_GeV",                         &fPion_MomY_Col_GeV,                         "fPion_MomY_Col_GeV/D");
  t1->Branch("Pion_MomZ_Col_GeV",                         &fPion_MomZ_Col_GeV,                         "fPion_MomZ_Col_GeV/D");

  t1->Branch("Neutron_Theta_Col",                         &fNeutron_Theta_Col,                         "fNeutron_Theta_Col/D");
  t1->Branch("Neutron_Phi_Col",                           &fNeutron_Phi_Col,                           "fNeutron_Phi_Col/D");
  t1->Branch("Neutron_Energy_Col_GeV",                    &fNeutron_Energy_Col_GeV,                    "fNeutron_Energy_Col_GeV/D");
  t1->Branch("Neutron_Mom_Col_GeV",                       &fNeutron_Mom_Col_GeV,                       "fNeutron_Mom_Col_GeV/D");
  t1->Branch("Neutron_MomX_Col_GeV",                      &fNeutron_MomX_Col_GeV,                      "fNeutron_MomX_Col_GeV/D");
  t1->Branch("Neutron_MomY_Col_GeV",                      &fNeutron_MomY_Col_GeV,                      "fNeutron_MomY_Col_GeV/D");
  t1->Branch("Neutron_MomZ_Col_GeV",                      &fNeutron_MomZ_Col_GeV,                      "fNeutron_MomZ_Col_GeV/D");

  // 6 Particles, electron, proton, scat. electron, photon, pion and neutron in proton's rest frame
  //-----------------------------------------------------------------------------------------------

  /* t1->Branch("Proton_Energy_RF_GeV",                     &fProton_Energy_RF_GeV,                     "fProton_Energy_RF_GeV/D"); */
  /* t1->Branch("Proton_Mom_RF_GeV",                        &fProton_Mom_RF_GeV,                        "fProton_Mom_RF_GeV/D"); */
  /* t1->Branch("Proton_MomZ_RF_GeV",                       &fProton_MomZ_RF_GeV,                       "fProton_MomZ_RF_GeV/D"); */
  /* t1->Branch("Proton_MomX_RF_GeV",                       &fProton_MomX_RF_GeV,                       "fProton_MomX_RF_GeV/D"); */
  /* t1->Branch("Proton_MomY_RF_GeV",                       &fProton_MomY_RF_GeV,                       "fProton_MomY_RF_GeV/D"); */

  /* t1->Branch("Electron_Theta_RF",                        &fElectron_Theta_RF,                        "fElectron_Theta_RF/D"); */
  /* t1->Branch("Electron_Phi_RF",                          &fElectron_Phi_RF,                          "fElectron_Phi_RF/D"); */
  /* t1->Branch("Electron_Energy_RF_GeV",                   &fElectron_Energy_RF_GeV,                   "fElectron_Energy_RF_GeV/D"); */
  /* t1->Branch("Electron_Mom_RF_GeV",                      &fElectron_Mom_RF_GeV,                      "fElectron_Mom_RF_GeV/D"); */
  /* t1->Branch("Electron_MomX_RF_GeV",                     &fElectron_MomX_RF_GeV,                     "fElectron_MomX_RF_GeV/D"); */
  /* t1->Branch("Electron_MomY_RF_GeV",                     &fElectron_MomY_RF_GeV,                     "fElectron_MomY_RF_GeV/D"); */
  /* t1->Branch("Electron_MomZ_RF_GeV",                     &fElectron_MomZ_RF_GeV,                     "fElectron_MomZ_RF_GeV/D"); */

  /* t1->Branch("ScatElec_Theta_RF",                        &fScatElec_Theta_RF,                        "fScatElec_Theta_RF/D"); */
  /* t1->Branch("ScatElec_Phi_RF",                          &fScatElec_Phi_RF,                          "fScatElec_Phi_RF/D"); */
  /* t1->Branch("ScatElec_Energy_RF_GeV",                   &fScatElec_Energy_RF_GeV,                   "fScatElec_Energy_RF_GeV/D"); */
  /* t1->Branch("ScatElec_Mom_RF_GeV",                      &fScatElec_Mom_RF_GeV,                      "fScatElec_Mom_RF_GeV/D"); */
  /* t1->Branch("ScatElec_MomX_RF_GeV",                     &fScatElec_MomX_RF_GeV,                     "fScatElec_MomX_RF_GeV/D"); */
  /* t1->Branch("ScatElec_MomY_RF_GeV",                     &fScatElec_MomY_RF_GeV,                     "fScatElec_MomY_RF_GeV/D"); */
  /* t1->Branch("ScatElec_MomZ_RF_GeV",                     &fScatElec_MomZ_RF_GeV,                     "fScatElec_MomZ_RF_GeV/D"); */

  /* t1->Branch("Photon_Energy_RF",                     &fPhoton_Energy_RF,                     "fPhoton_Energy_RF/D"); */
  /* t1->Branch("Photon_Mom_RF",                        &fPhoton_Mom_RF,                        "fPhoton_Mom_RF/D"); */
  /* t1->Branch("Photon_MomX_RF",                       &fPhoton_MomX_RF,                       "fPhoton_MomX_RF/D"); */
  /* t1->Branch("Photon_MomY_RF",                       &fPhoton_MomY_RF,                       "fPhoton_MomY_RF/D"); */
  /* t1->Branch("Photon_MomZ_RF",                       &fPhoton_MomZ_RF,                       "fPhoton_MomZ_RF/D"); */

  /* t1->Branch("Photon_Theta_RF",                          &fPhoton_Theta_RF,                          "fPhoton_Theta_RF/D"); */
  /* t1->Branch("Photon_Phi_RF",                            &fPhoton_Phi_RF,                            "fPhoton_Phi_RF/D"); */
  /* t1->Branch("Photon_Energy_RF_GeV",                     &fPhoton_Energy_RF_GeV,                     "fPhoton_Energy_RF_GeV/D"); */
  /* t1->Branch("Photon_Mom_RF_GeV",                        &fPhoton_Mom_RF_GeV,                        "fPhoton_Mom_RF_GeV/D"); */
  /* t1->Branch("Photon_MomX_RF_GeV",                       &fPhoton_MomX_RF_GeV,                       "fPhoton_MomX_RF_GeV/D"); */
  /* t1->Branch("Photon_MomY_RF_GeV",                       &fPhoton_MomY_RF_GeV,                       "fPhoton_MomY_RF_GeV/D"); */
  /* t1->Branch("Photon_MomZ_RF_GeV",                       &fPhoton_MomZ_RF_GeV,                       "fPhoton_MomZ_RF_GeV/D"); */

  /* t1->Branch("Pion_Theta_RF",                            &fPion_Theta_RF,                            "fPion_Theta_RF/D"); */
  /* t1->Branch("Pion_Phi_RF",                              &fPion_Phi_RF,                              "fPion_Phi_RF/D"); */
  /* t1->Branch("Pion_Energy_RF_GeV",                       &fPion_Energy_RF_GeV,                       "fPion_Energy_RF_GeV/D"); */
  /* t1->Branch("Pion_Mom_RF_GeV",                          &fPion_Mom_RF_GeV,                          "fPion_Mom_RF_GeV/D"); */
  /* t1->Branch("Pion_MomX_RF_GeV",                         &fPion_MomX_RF_GeV,                         "fPion_MomX_RF_GeV/D"); */
  /* t1->Branch("Pion_MomY_RF_GeV",                         &fPion_MomY_RF_GeV,                         "fPion_MomY_RF_GeV/D"); */
  /* t1->Branch("Pion_MomZ_RF_GeV",                         &fPion_MomZ_RF_GeV,                         "fPion_MomZ_RF_GeV/D"); */

  /* t1->Branch("Neutron_Theta_RF",                         &fNeutron_Theta_RF,                         "fNeutron_Theta_RF/D"); */
  /* t1->Branch("Neutron_Phi_RF",                           &fNeutron_Phi_RF,                           "fNeutron_Phi_RF/D"); */
  /* t1->Branch("Neutron_Energy_RF_GeV",                    &fNeutron_Energy_RF_GeV,                    "fNeutron_Energy_RF_GeV/D"); */
  /* t1->Branch("Neutron_Mom_RF_GeV",                       &fNeutron_Mom_RF_GeV,                       "fNeutron_Mom_RF_GeV/D"); */
  /* t1->Branch("Neutron_MomX_RF_GeV",                      &fNeutron_MomX_RF_GeV,                      "fNeutron_MomX_RF_GeV/D"); */
  /* t1->Branch("Neutron_MomY_RF_GeV",                      &fNeutron_MomY_RF_GeV,                      "fNeutron_MomY_RF_GeV/D"); */
  /* t1->Branch("Neutron_MomZ_RF_GeV",                      &fNeutron_MomZ_RF_GeV,                      "fNeutron_MomZ_RF_GeV/D"); */

  // ----------------------------------------------------------------------------------------------------------------------------------

  /* t1->Branch("Photon_Theta_RF",                           &fPhoton_Theta_RF,                           "fPhoton_Theta_RF/D"); */
  /* t1->Branch("Photon_Theta_Col",                          &fPhoton_Theta_Col,                          "fPhoton_Theta_Col/D"); */
  /* t1->Branch("NRecorded",                                 &fNRecorded,                                 "fNRecorded/I"); */
  /* t1->Branch("NGenerated",                                &fNGenerated,                                "fNGeneratedd/I"); */

  t1->Branch("Epsilon",                                   &fEpsilon,                                   "fEpsilon/D");
  t1->Branch("Phi",                                       &fPhi,                                       "fPhi/D");
  t1->Branch("PhiS",                                      &fPhiS,                                      "fPhiS/D");
  t1->Branch("W_GeV",                                     &fW_GeV,                                     "fW_GeV/D");
  t1->Branch("W_Prime_GeV",                               &fW_Prime_GeV,                               "fW_Prime_GeV/D");

  t1->Branch("Qsq_GeV",                                   &fQsq_GeV,                                   "fQsq_GeV/D");
  t1->Branch("T_Para_GeV",                                &fT_Para_GeV,                                "fT_Para_GeV/D");
  t1->Branch("T_GeV",                                     &fT_GeV,                                     "fT_GeV/D");
  t1->Branch("x",                                         &fx,                                         "fx/D");
  t1->Branch("y",                                         &fy,                                         "fy/D");
  t1->Branch("z",                                         &fz,                                         "fz/D");

  t1->Branch("Flux_Factor_RF",                            &fFlux_Factor_RF,                            "fFlux_Factor_RF/D");
  t1->Branch("Flux_Factor_Col",                           &fFlux_Factor_Col,                           "fFlux_Factor_Col/D");
  t1->Branch("Jacobian_CM",                               &fJacobian_CM,                               "fJacobian_CM/D");
  t1->Branch("Jacobian_CM_RF",                            &fJacobian_CM_RF,                            "fJacobian_CM_RF/D");
  t1->Branch("Jacobian_CM_Col",                           &fJacobian_CM_Col,                           "fJacobian_CM_Col/D");
  t1->Branch("EventWeight",                               &fEventWeight,                               "fEventWeight/D");
  t1->Branch("Sigma_Col",                                 &fSigma_Col,                                 "fSigma_Col/D");
  t1->Branch("Sig_VR"  ,                                  &fSig_VR,                                    "fSig_VR/D");
  t1->Branch("Sig_L",                                     &fSig_L,                                     "fSig_L/D");
  t1->Branch("Sig_T",                                     &fSig_T,                                     "fSig_T/D");

  t1->Branch("A",                                         &fA,                                         "fA/D");
  t1->Branch("Vertex_X",                                  &fVertex_X,                                  "fVertex_X/D");
  t1->Branch("Vertex_Y",                                  &fVertex_Y,                                  "fVertex_Y/D");
  t1->Branch("Vertex_Z",                                  &fVertex_Z,                                  "fVertex_Z/D");
  
  /* t1->Branch("Pion_Energy_CM_GeV",                        &fPion_Energy_CM_GeV,                        "fPion_Energy_CM_GeV/D"); */
  /* t1->Branch("Pion_Mom_CM_GeV",                           &fPion_Mom_CM_GeV,                           "fPion_Mom_CM_GeV/D"); */
  /* t1->Branch("BetaX_Col_RF",                              &fBetaX_Col_RF,                              "fBetaX_Col_RF/D"); */
  /* t1->Branch("BetaY_Col_RF",                              &fBetaY_Col_RF,                              "fBetaY_Col_RF/D"); */
  /* t1->Branch("BetaZ_Col_RF",                              &fBetaZ_Col_RF,                              "fBetaZ_Col_RF/D"); */
  /* t1->Branch("Beta_Col_RF",                               &fBeta_Col_RF,                               "fBeta_Col_RF/D"); */
  /* t1->Branch("Gamma_Col_RF",                              &fGamma_Col_RF,                              "fGamma_Col_RF/D"); */
  /* t1->Branch("Beta_CM_RF",                                &fBeta_CM_RF,                                "fBeta_CM_RF/D"); */
  /* t1->Branch("Gamma_CM_RF",                               &fGamma_CM_RF,                               "fGamma_CM_RF/D"); */

  /* t1->Branch("XMomConserve",                              &fXMomConserve,                              "fXMomConserve/D"); */
  /* t1->Branch("YMomConserve",                              &fYMomConserve,                              "fYMomConserve/D"); */
  /* t1->Branch("ZMomConserve",                              &fZMomConserve,                              "fZMomConserve/D"); */
  /* t1->Branch("EnergyConserve",                            &fEnergyConserve,                            "fEnergyConserve/D"); */

  /* t1->Branch("XMomConserve_RF",                           &fXMomConserve_RF,                           "fXMomConserve_RF/D"); */
  /* t1->Branch("YMomConserve_RF",                           &fYMomConserve_RF,                           "fYMomConserve_RF/D"); */
  /* t1->Branch("ZMomConserve_RF",                           &fZMomConserve_RF,                           "fZMomConserve_RF/D"); */
  /* t1->Branch("EnergyConserve_RF",                         &fEnergyConserve_RF,                         "fEnergyConserve_RF/D"); */
  /* t1->Branch("testsig",                                   &ftestsig,                                   "ftestsig/D"); */

}
