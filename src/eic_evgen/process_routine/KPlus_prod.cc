#include "reaction_routine.h"
#include "eic.h"

using namespace std;

KPlus_Production::KPlus_Production() { 

  cout << "Program Start" << endl;

}

/// KPlus_Production 
// 09/02/22 - SJDK - Added new hadron string argument
KPlus_Production::KPlus_Production(TString particle_str, TString hadron_str) { 

  rParticle = particle_str;
  rHadron = hadron_str;

}

KPlus_Production::~KPlus_Production() {

  ppiOut.close();
  ppiDetails.close();

}

void KPlus_Production::process_reaction() {
 
  Init();

  if (gOutputType == "Pythia6"){
    KPlus_Pythia6_Out_Init();
  }
  else if (gOutputType == "HEPMC3"){
    KPlus_HEPMC3_Out_Init();
  }

  for( long long int i = 0; i < rNEvents; i++ ) {
 
    rNEvent_itt = i;
    fNGenerated ++;
 
    Progress_Report();  // This is happens at each 10% of the total event is processed
    Processing_Event();
  }
 
  Detail_Output();
 
}

void KPlus_Production::Init() {
  
  pim* myPim;
  
  pd = dynamic_cast<pim*>(myPim);
  rParticle_charge = ExtractCharge(rParticle);

  sTFile = Form("./LundFiles/eic_%s.txt", gfile_name.Data());
  sLFile= Form("./LundFiles/eic_input_%s.dat", gfile_name.Data());
   
  ppiOut.open( sLFile.c_str() );
  ppiDetails.open( sTFile.c_str() );
	
  qsq_ev = 0, t_ev = 0, w_neg_ev = 0, w_ev = 0;
  rNEvents = fNEvents;
  rNEvent_itt = 0;

  // 02/06/21 SJDK 
  // Set these values once the beam energies are read in
  fPSF = ( fEBeam * ( fScatElec_E_Hi - fScatElec_E_Lo ) *( sin( fScatElec_Theta_F ) - sin( fScatElec_Theta_I ) ) * 2 * fPI *( sin( fPion_Theta_F     ) - sin( fPion_Theta_I     ) ) * 2 * fPI );
  fElectron_Kin_Col_GeV = fEBeam;
  fElectron_Kin_Col = fElectron_Kin_Col_GeV * 1000.0;

  // cout << rNEvents << "    " << fNEvents << endl;
	
  rFermiMomentum = pd->fermiMomentum();

  // Proton in collider (lab) frame

  r_lproton = GetProtonVector_lab();
  r_lprotong = GetProtonVector_lab() * fm;

  /// Electron in collider (lab) frame

  cout << "Fermi momentum: " << rFermiMomentum << endl;

  r_lelectron	 = GetElectronVector_lab();
  r_lelectrong = r_lelectron * fm;

  ///*--------------------------------------------------*/
  /// Getting the particle mass from the data base
  /// Particle X is the produced meson, in this case, the K+
  ///*--------------------------------------------------*/
  produced_X = ParticleEnum(rParticle);
  fX_Mass = ParticleMass(produced_X)*1000; //MeV
  fX_Mass_GeV = fX_Mass/1000; //GeV

  cout << rParticle << "  " << produced_X << "  " << fX_Mass_GeV <<  endl;
  cout << rParticle_charge << endl;

  if(rHadron == "Lambda"){
    rParticle_scat_hadron  = "Lambda"; 
    recoil_hadron  = Lambda;
    f_Scat_hadron_Mass = fLambda_Mass;
    f_Scat_hadron_Mass_GeV = f_Scat_hadron_Mass/1000;
    cout<<"Particle = "<<rParticle_scat_hadron<<" with mass = "<< f_Scat_hadron_Mass << endl;
  }
  else if (rHadron == "Sigma0"){
    rParticle_scat_hadron  = "Sigma0";
    recoil_hadron  = Sigma0;
    f_Scat_hadron_Mass  = fSigma_Mass;
    f_Scat_hadron_Mass_GeV = f_Scat_hadron_Mass/1000;
    cout<<"Particle = "<<rParticle_scat_hadron<<" with mass = "<< f_Scat_hadron_Mass << endl;
  }
  else{
    cout << "Error, expect Lambda or Sigma0 specified as hadron for K+ production module." << endl;
    cout << "Defaulting to Lambda production." << endl;
    rParticle_scat_hadron  = "Lambda"; 
    recoil_hadron  = Lambda;
    f_Scat_hadron_Mass = fLambda_Mass;
    f_Scat_hadron_Mass_GeV = f_Scat_hadron_Mass/1000;
    cout<<"Particle = "<<rParticle_scat_hadron<<" with mass = "<< f_Scat_hadron_Mass << endl;
  }

  rDEG2RAD   = fPI/180.0;

  fX_Theta_I = 0.0 * rDEG2RAD ;
  fX_Theta_F = 43.0 * rDEG2RAD;

  cout << "Produced particle in exclusive production: " << rParticle << ";  with mass: " << fX_Mass << " MeV "<< endl;
  cout << fEBeam << " GeV electrons on " << fPBeam << " GeV ions" << endl;

  // Depending upon beam energy combination, set the value for the max weight from the non normalised version to then generate unit weights
  // The values were determined from a set of 100 x 1B events thrown runs, the mean weight value + 6.5 sigma was taken as the "max" weight for a given beam energy combination
  // Probability of being more than 6.5 sigma away is over 1 in 12.5B
  // The weight has to be scaled by the number thrown in the current calculation 
  // fEventWeight is now independent of the number of events thrown

  // SJDK 21/06/21 - Commented out for now, reverting to old weighting method

  // if ((fEBeam == 5.0 ) && (fPBeam == 41.0) ){
  //   //fEventWeightCeil = 0.0221836 * (1000000000); // Old value
  //   fEventWeightCeil = 0.002296 * (1000000000);
  // }

  // else if ((fEBeam == 5.0 ) && (fPBeam == 100.0) ){
  //   //fEventWeightCeil = 0.30281 * (1000000000); // Old value
  //   fEventWeightCeil = 0.023960 * (1000000000);
  // }

  // else if ((fEBeam == 10.0 ) && (fPBeam == 100.0) ){
  //   //fEventWeightCeil = 1.77775 * (1000000000); // Old value
  //   fEventWeightCeil = 0.201569 * (1000000000);
  // }
  // else {
  //   fEventWeightCeil = 1.0 * (100000000);
  //   cout << endl << "!!!!! WARNING !!!!!" << endl;
  //   cout << "Beam energy combination not recognised, weight ceiling set to 1." << endl;
  //   cout << "!!!!! WARNING !!!!!" << endl << endl;
  // }

}

/*--------------------------------------------------*/

void KPlus_Production::Processing_Event() {

  // ----------------------------------------------------
  // Considering Fermi momentum for the proton
  // ----------------------------------------------------

  if( kCalcFermi ) {
    Consider_Proton_Fermi_Momentum();
 
  }

  // ----------------------------------------------------
  // Boost vector from collider (lab) frame to protons rest frame (Fix target)
  // ----------------------------------------------------
 
  beta_col_rf = r_lproton.BoostVector();        
  fGamma_Col_RF = 1.0/sqrt( 1 - pow( beta_col_rf.Mag() , 2 ) );

  // ---------------------------------------------------------------------
  // Specify the energy and solid angle of scatterd electron in Collider (lab) frame
  // ---------------------------------------------------------------------
  fScatElec_Theta_Col  = acos( fRandom->Uniform( cos( fScatElec_Theta_I ) , cos( fScatElec_Theta_F ) ) );
  fScatElec_Phi_Col    = fRandom->Uniform( 0 , 2.0 * fPi);
  fScatElec_Energy_Col = fRandom->Uniform( fScatElec_E_Lo * fElectron_Energy_Col , fScatElec_E_Hi * fElectron_Energy_Col );

  // ----------------------------------------------------
  // Produced Particle X in Collider frame
  // ----------------------------------------------------  


  /// The generic produced particle in the exclusive reaction is labelled as X 
  fX_Theta_Col      = acos( fRandom->Uniform( cos(fX_Theta_I), cos(fX_Theta_F ) ) ); 
  fX_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi );
    	
  // ---------------------------------------------------------------------
  // Specify the energy and solid angle of scatterd electron in Collider (lab) frame
  // ---------------------------------------------------------------------

  fScatElec_Mom_Col  = sqrt( pow( fScatElec_Energy_Col,2) - pow( fElectron_Mass , 2) );
  fScatElec_MomZ_Col = ( fScatElec_Mom_Col * cos(fScatElec_Theta_Col) );  
  fScatElec_MomX_Col = ( fScatElec_Mom_Col * sin(fScatElec_Theta_Col) * cos(fScatElec_Phi_Col) );
  fScatElec_MomY_Col = ( fScatElec_Mom_Col * sin(fScatElec_Theta_Col) * sin(fScatElec_Phi_Col) );

  r_lscatelec.SetPxPyPzE( fScatElec_MomX_Col, fScatElec_MomY_Col, fScatElec_MomZ_Col, fScatElec_Energy_Col);
 
  r_lscatelecg = r_lscatelec * fm;

  // ----------------------------------------------------
  // Photon in collider (lab) frame and Qsq
  // ----------------------------------------------------

  r_lphoton  = r_lelectron - r_lscatelec;
  r_lphotong = r_lelectrong - r_lscatelecg;

  fQsq_GeV = -1.* r_lphotong.Mag2();

  if ( fQsq_GeV < 3.0 ) {
    qsq_ev++;
    return;
  }

  // ----------------------------------------------------
  // W square, Invariant Mass (P_g + P_p)^2
  // ----------------------------------------------------
 
  TLorentzVector lwg;
  lwg = r_lprotong + r_lphotong;
  fW_GeV    = lwg.Mag();
  fWSq_GeV  = lwg.Mag2();
    
  if ( fWSq_GeV < 0 ) { 
    w_neg_ev++;
    return;
  }    
 
  // ---------------------------------------------------------
  // Kaon momentum in collider frame, analytic solution starts
  // ---------------------------------------------------------
 
  double fupx = sin( fX_Theta_Col ) * cos( fX_Phi_Col );
  double fupy = sin( fX_Theta_Col ) * sin( fX_Phi_Col );
  double fupz = cos( fX_Theta_Col );
 
  double fuqx = sin( r_lphoton.Theta() ) * cos( r_lphoton.Phi() );
  double fuqy = sin( r_lphoton.Theta() ) * sin( r_lphoton.Phi() );
  double fuqz = cos( r_lphoton.Theta() );
 
  double fa = -(r_lphoton.Vect()).Mag() * ( fupx * fuqx +  fupy * fuqy +  fupz * fuqz );
  double fb = pow ( (r_lphoton.Vect()).Mag() , 2 );
  double fc = r_lphoton.E() + fProton_Mass;
 
  fa = ( fa - std::abs( (r_lproton.Vect()).Mag() ) * ( ( ( r_lproton.X() / (r_lproton.Vect()).Mag() ) * fupx ) + 
 						       ( ( r_lproton.Y() / (r_lproton.Vect()).Mag() ) * fupy ) + 
 						       ( ( r_lproton.Z() / (r_lproton.Vect()).Mag() ) * fupz ) ) );
     
  double factor = ( pow( (r_lproton.Vect()).Mag() , 2 ) + 2.0 * (r_lphoton.Vect()).Mag() * (r_lproton.Vect()).Mag() *  
		    ( ( ( r_lproton.X() / (r_lproton.Vect()).Mag() ) * fuqx ) + 
		      ( ( r_lproton.Y() / (r_lproton.Vect()).Mag() ) * fuqy ) + 
		      ( ( r_lproton.Z() / (r_lproton.Vect()).Mag() ) * fuqz ) ) );
     
  fb =  fb + factor;  
  fc = r_lphoton.E() + r_lproton.E();
     
  double ft = fc * fc - fb + fX_Mass * fX_Mass - fProton_Mass * fProton_Mass;
    
  double fQA = 4.0 * ( fa * fa - fc * fc );
  double fQB = 4.0 * fc * ft;

  double fQC = -4.0 * fa * fa * fX_Mass * fX_Mass - ft * ft;    
 
  fradical = fQB * fQB - 4.0 * fQA * fQC;
 
  fepi1 = ( -fQB - sqrt( fradical ) ) / ( 2.0 * fQA );
  fepi2 = ( -fQB + sqrt( fradical ) ) / ( 2.0 * fQA );
 
  // ---------------------------------------------------------
  // Particle X momentum in collider frame, analytic solution ends
  // ---------------------------------------------------------
         
  r_lX.SetPxPyPzE( (sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * sin(fX_Theta_Col) * cos(fX_Phi_Col),
		   ( sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * sin(fX_Theta_Col) * sin(fX_Phi_Col),
		   ( sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * cos(fX_Theta_Col),
		   fepi1 );

  r_lX_g = r_lX * fm;

  // ----------------------------------------------------
  // Scattered proton collider (lab) frame
  // ----------------------------------------------------
  
  r_l_scat_hadron.SetPxPyPzE( ( r_lproton + r_lelectron - r_lscatelec - r_lX).X(),
			       ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Y(),
			       ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Z(),
 			       sqrt( pow( ( ( ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Vect() ).Mag()),2) +
  				     pow( f_Scat_hadron_Mass ,2 ) ) );

  r_l_scat_hadron_g = r_l_scat_hadron * fm;

  // ----------------------------------------------------------------------------------------------
  // Calculate w = (proton + photon)^2
  // ----------------------------------------------------------------------------------------------
     
  if ( fW_GeV < 3.0 || fW_GeV > 10.6 ) {
    w_ev++;
    return;
  }

  r_lw = r_lproton + r_lphoton;
  fW = r_lw.Mag();


  // ----------------------------------------------------------------------------------------------
  // Calculate w prime w' = (proton + photon - pion)^2                                             
  // ----------------------------------------------------------------------------------------------

  lwp = r_lprotong + r_lphotong - r_lX_g;
  fW_Prime_GeV = lwp.Mag();    

  fsini = r_lelectron + r_lproton;
  fsfin = r_lscatelec + r_lX + r_l_scat_hadron;
     
  fsinig = fsini * fm;
  fsfing = fsfin * fm; 
  // SJDK 15/06/21 - Mandlestam S conservation check - doesn't actually seem to be utilised?
  fMandSConserve = std::abs( fsinig.Mag() - fsfing.Mag() );

  // SJDK 15/06/21 - Added integer counters for conservation law check and for NaN check
  if (r_lX.E() != r_lX.E()){ // SJDK 15/06/21 - If the energy of the produced meson is not a number, return and add to counter
    fNaN++;
    return;
  }  
  kSConserve = false;
  if( std::abs( fsinig.Mag() - fsfing.Mag() ) < fDiff ) {
    kSConserve = true;
  }
  // For the Kaon case, the energy difference threshold is set higher for now
  if ( pd->CheckLaws( r_lelectron, r_lproton, r_lscatelec, r_lX, r_l_scat_hadron, 10) != 1 ){
    fConserve++;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  //                                          Start                                         //
  // Transformation of e', pi- and recoil proton to target's rest frame without energy loss //
  ////////////////////////////////////////////////////////////////////////////////////////////
 
  lproton_rf = r_lproton;
  lproton_rf.Boost(-beta_col_rf);
  lproton_rfg = lproton_rf * fm;
 
  lelectron_rf = r_lelectron;
  lelectron_rf.Boost(-beta_col_rf);
  lelectron_rfg = lelectron_rf * fm;
 
  lscatelec_rf = r_lscatelec;
  lscatelec_rf.Boost(-beta_col_rf);
  lscatelec_rfg = lscatelec_rf * fm;
     
  lphoton_rf = r_lphoton;
  lphoton_rf.Boost(-beta_col_rf);
  lphoton_rfg = lphoton_rf * fm;
     
  lX_rf = r_lX;
  lX_rf.Boost(-beta_col_rf);
  lX_rfg = lX_rf * fm;
         
  l_scat_hadron_rf = r_l_scat_hadron;
  l_scat_hadron_rf.Boost(-beta_col_rf);
  l_scat_hadron_rf_g = l_scat_hadron_rf * fm;

  ////////////////////////////////////////////////////////////////////////////////////////////
  //                                          End                                           //
  // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
  ////////////////////////////////////////////////////////////////////////////////////////////


  // -----------------------------------------------------------------------------------------
  // Calculate -t
  // -----------------------------------------------------------------------------------------

  fBeta_CM_RF           = (lphoton_rf.Vect()).Mag() / ( lphoton_rf.E() + fProton_Mass );
  fGamma_CM_RF          = ( lphoton_rf.E() + fProton_Mass ) / fW;
  fX_Energy_CM       = ( pow( fW , 2) + pow(fX_Mass , 2) - pow(f_Scat_hadron_Mass , 2) ) / ( 2.0 * fW);    
  fX_Mom_CM          = sqrt( pow(fX_Energy_CM , 2) - pow(fX_Mass , 2));    
  fX_Energy_CM_GeV   = fX_Energy_CM / 1000.0;
  fX_Mom_CM_GeV      = fX_Mom_CM / 1000.0;

  // this equation is valid for parallel kinematics only!
  fT_Para = ( pow(((r_lphoton.Vect()).Mag() - (r_lX.Vect()).Mag()),2) - pow((r_lphoton.E() - r_lX.E()),2));
  fT_Para_GeV = fT_Para/1000000.0;

  lt = r_lphoton - r_lX;
  ltg = lt * fm;
 
  fT = -1.*lt.Mag2();
  fT_GeV = -1.*ltg.Mag2();

  if ( gKinematics_type == 1 && fT_GeV > 1.2 ) {
    t_ev++;
    return;
  }
     
  if ( gKinematics_type == 2 && fT_GeV > 1.3 ) {
    t_ev++;
    return;
  }

  fx = fQsq_GeV / ( 2.0 * r_lprotong.Dot( r_lphotong ) );
  fy = r_lprotong.Dot( r_lphotong ) / r_lprotong.Dot( r_lelectrong );
  fz = r_lX.E()/r_lphoton.E();    

  // -------------------------------------------------------------------------------------------------------
  // Calculation of Phi  ( azimuthal angle of pion momentum w.r.t lepton plane in target's rest frame)
  // Calculation of PhiS ( azimuthal angle of target polarization w.r.t lepton plane in target's rest frame)
  // -------------------------------------------------------------------------------------------------------  

  v3Photon.SetX( lphoton_rfg.X() );     
  v3Photon.SetY( lphoton_rfg.Y() );     
  v3Photon.SetZ( lphoton_rfg.Z() );    

  v3Electron.SetX( lelectron_rfg.X() ); 
  v3Electron.SetY( lelectron_rfg.Y() ); 
  v3Electron.SetZ( lelectron_rfg.Z() );

  v3X.SetX( lX_rfg.X() ) ;        
  v3X.SetY( lX_rfg.Y() ) ;        
  v3X.SetZ( lX_rfg.Z() );

  v3S.SetX( -1 );                       
  v3S.SetY( 0 );                        
  v3S.SetZ( 0 );        

  v3PhotonUnit = v3Photon.Unit();    
  v3QxL        = v3Photon.Cross(v3Electron);
  v3QxP        = v3Photon.Cross(v3X);
  v3QxS        = v3Photon.Cross(v3S);
  v3LxP        = v3Electron.Cross(v3X);
  v3LxS        = v3Electron.Cross(v3S);
  v3PxL        = v3X.Cross(v3Electron);
  v3QUnitxL    = v3PhotonUnit.Cross(v3Electron);
  v3QUnitxP    = v3PhotonUnit.Cross(v3X);
  v3QUnitxS    = v3PhotonUnit.Cross(v3S);

  // Get the Phi scattering angle with respect to the electron scattering plane
  fPhi   = Get_Phi_X_LeptonPlane_RF ();

  // Get the Phi scattering angle with respect to the electron scattering plane
  fPhiS  = Get_Phi_TargPol_LeptonPlane_RF();

  fTheta_X_Photon_RF       = fRAD2DEG * acos( ( v3Photon.Dot( v3X     ) ) / ( v3Photon.Mag()  * v3X.Mag()    ) );
  if ( fTheta_X_Photon_RF < 0 ) { fTheta_X_Photon_RF = 180.0 + fTheta_X_Photon_RF; }

  // -----------------------------------------------------------------------------------
  // If we have fermi momentum then epsilon should be in rest frame 
  // The theta angle of scattered angle used in expression of epsilon is the angle 
  // with respect to direction of incoming electron in the rest frame of target nucleon
  // epsilon=1./(1.+ 2.*(pgam_restg**2)/q2g * *(tand(thscat_rest/2.))**2)
  // -----------------------------------------------------------------------------------
 
  double fTheta_EEp = (lelectron_rf.Vect()).Angle(lscatelec_rf.Vect());

  fEpsilon = 1.0 / ( 1.0 + 2.0 * ( pow( (lphoton_rfg.Vect()).Mag(),2)/fQsq_GeV ) * pow( tan( fTheta_EEp / 2 ) , 2 ) );

  // ----------------------------------------------------
  // Virtual Photon flux factor in units of 1/(GeV*Sr)
  // ----------------------------------------------------
  fFlux_Factor_Col = (fAlpha/(2.0*pow(fPi,2))) * (r_lscatelecg.E() / r_lelectrong.E()) * 
    ( pow(fW_GeV,2) - pow(fProton_Mass_GeV,2) ) / (2.0*fProton_Mass_GeV*fQsq_GeV*(1.0 - fEpsilon));
         
  fFlux_Factor_RF = ( fAlpha / ( 2.0 * pow( fPi , 2 ) ) ) * ( lscatelec_rfg.E() / lelectron_rfg.E() ) *
    ( pow( fW_GeV , 2 ) - pow( fProton_Mass_GeV , 2 ) ) /
    ( 2.0 * fProton_Mass_GeV * fQsq_GeV * ( 1.0 - fEpsilon ) );

  // ----------------------------------------------------
  //  Jacobian  dt/dcos(theta*)dphi in units of GeV2/sr
  // ----------------------------------------------------
  fJacobian_CM = ( (lphoton_rfg.Vect()).Mag() - fBeta_CM_RF * lphoton_rfg.E() ) / ( fGamma_CM_RF * ( 1.0 - pow(fBeta_CM_RF,2) ) );
 
  fA = fJacobian_CM * fX_Mom_CM_GeV / fPi;
 
  // ----------------------------------------------------
  // Jacobian dOmega* / dOmega dimensionless
  // ----------------------------------------------------
  fJacobian_CM_RF  = ( pow((lX_rf.Vect()).Mag(),2)*fW) / 
    ( fX_Mom_CM * std::abs( ( fProton_Mass + lphoton_rf.E()) * (lX_rf.Vect()).Mag() - 
			    ( lX_rf.E() * (lphoton_rf.Vect()).Mag() * cos( lX_rf.Theta() ) ) ) );
 
  fJacobian_CM_Col = ( ( pow((r_lX.Vect()).Mag(),2) * fW ) /
		       ( fX_Mom_CM * std::abs( ( fProton_Mass + r_lphoton.E() ) * (r_lX.Vect()).Mag() -
					       ( r_lX.E() * (r_lphoton.Vect()).Mag() * cos( r_lX.Theta() ) ) ) ) );


  //	 cout <<  lX_rf.Vect().Mag() << "  " << << << << << << << << endl;
  //	 cout << fJacobian_CM_RF << "    " << fJacobian_CM_Col << endl;

  // -----------------------------------------------------------------------------------------------------------
  // CKY sigma L and T starts
  // -----------------------------------------------------------------------------------------------------------

  //	 r_fSig_T = 1;
  //	 r_fSig_L = 1;
  // 
  // -------------------------------------------------------------------------------------------
 
  r_fSig = Get_Total_Cross_Section();
 
  // -----------------------------------------------------------------------------------------------------------
  // CKY sigma L and T ends
  // -----------------------------------------------------------------------------------------------------------
 
  fSigma_Col = r_fSig * fFlux_Factor_Col * fA * fJacobian_CM_Col;

  if ( ( fSigma_Col <= 0 ) || std::isnan( fSigma_Col ) ) { 
    fNSigmaNeg ++;
    return;
  }
     
  // -----------------------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------------------
  //             Lab cross section     Phase Space   Conversion     Luminosity                Total events tried
  // Hz        = ub / ( sr^2 * GeV ) * GeV * sr^2 * ( cm^2 / ub ) * ( # / ( cm^2 * sec ) ) / ( # )

  // SJDK 11/05/21 - This is the previous non unit weight
  // SJDK 24/06/21 - Explicitly taking the absolute value of the weigth such that the value is positive!
  fEventWeight = abs(fSigma_Col * fPSF * fuBcm2 * fLumi / fNEvents);   // in Hz
  // SJDK 21/06/21 - Commenting out "unit weight" calculation for now, reverting to old version
  // SJDK 11/05/21 - New weight calculation, division by ceiling weight value to attempt to get a "unit" value
  //fEventWeight = abs(fSigma_Col * fPSF * fuBcm2 * fLumi )/fEventWeightCeil;
  // if ( (fEventWeight > 1) || (fEventWeight <= 0) ){
  //   fNWeightUnphys ++;
  //   return;
  // }
  // SJDK 21/06/21 - Reversed sign of condition here, actually want to reject those with a value less than the random number
  // fEventWeightRn = fRandom->Uniform( 0, 1.0);
  // if ( fEventWeight < fEventWeightRn ){
  //   fNWeightReject ++;
  //   return;
  // }
  fNRecorded++;
  fLundRecorded++;
  fRatio = fNRecorded / fNGenerated;

  if (gOutputType == "Pythia6"){
      KPlus_Pythia6_Output();
  }
  else if (gOutputType == "LUND"){
    Lund_Output();
  }
  else if (gOutputType == "HEPMC3"){
    KPlus_HEPMC3_Output();
  }	       
}

void KPlus_Production::Progress_Report() {

  dFractTime = time(0);

  if ( rNEvent_itt % ( rNEvents / 10 ) == 0 ) {
    cout << "Event: " << setw(8) << rNEvent_itt 
	 << "     % of events " << setw(4) << ((1.0*rNEvent_itt)/(1.0*rNEvents))*100.0
	 << "   Day: " <<  dFractTime.GetDay() 
	 << "   Time:   " << dFractTime.GetHour() 
	 << ":" << dFractTime.GetMinute() 
	 << ":" << dFractTime.GetSecond() 
	 << endl;	  
  }
}

TLorentzVector KPlus_Production::GetProtonVector_lab() {

  // Crossing angle
  //	 fProton_Theta_Col = 0.050;
  //	 fProton_Theta_Col = 0.025;
  // Set crossing angle to 0 for fun4all, also required for ATHENA simulations
  fProton_Theta_Col = 0.0;

  //     fProton_Phi_Col   = fPi; 
  fProton_Phi_Col   = fProton_incidence_phi; 

  fProton_Mom_Col   = fPBeam * 1e3; 
  fVertex_X         = 0.; 
  fVertex_Y         = 0.; 
  fVertex_Z         = 0.; 

 
  TLorentzVector lproton( fProton_Mom_Col * sin(fProton_Theta_Col) * cos(fProton_Phi_Col),
			  fProton_Mom_Col * sin(fProton_Theta_Col) * sin(fProton_Phi_Col),
			  fProton_Mom_Col * cos(fProton_Theta_Col),
			  sqrt( pow( fProton_Mom_Col , 2 ) + pow( fProton_Mass , 2 ) ) ); 

  return lproton;

}

//*--------------------------------------------------*/
// Proton in collider (lab) frame
// ----------------------------------------------------

void KPlus_Production::Consider_Proton_Fermi_Momentum() {

  fProton_Mom_Col   = fProton_Mom_Col + rFermiMomentum;
  fProton_Theta_Col = acos( fRandom->Uniform( cos(0.0) , cos(fPi) ) );
  fProton_Phi_Col   = fRandom->Uniform( 0 , 360 );

  double px, py, pz, e;

  px = fProton_Mom_Col * sin(fProton_Theta_Col) * cos(fProton_Phi_Col);
  py = fProton_Mom_Col * sin(fProton_Theta_Col) * sin(fProton_Phi_Col);
  pz = fProton_Mom_Col * cos(fProton_Theta_Col);
  e  = sqrt( pow( fProton_Mom_Col , 2 ) + pow( fProton_Mass , 2 ) );

  r_lproton.SetPxPyPzE(px,py,pz,e);

  r_lprotong = r_lproton*fm;

}

// ----------------------------------------------------
// Electron in collider (lab) frame
// ----------------------------------------------------

TLorentzVector KPlus_Production::GetElectronVector_lab() {

  fElectron_Energy_Col = fElectron_Kin_Col; 
  fElectron_Mom_Col    = sqrt( pow(fElectron_Energy_Col , 2) - pow(fElectron_Mass , 2) );
  fElectron_Theta_Col  = fPi;
  fElectron_Phi_Col    = 0.0;
  fElectron_MomZ_Col   = fElectron_Mom_Col * cos(fElectron_Theta_Col);  
  fElectron_MomX_Col   = fElectron_Mom_Col * sin(fElectron_Theta_Col) * cos(fElectron_Phi_Col);
  fElectron_MomY_Col   = fElectron_Mom_Col * sin(fElectron_Theta_Col) * sin(fElectron_Phi_Col);  

  cout << "Define: " << fElectron_MomZ_Col << "    "<< fElectron_Mom_Col << "  " << cos(fElectron_Theta_Col) << endl;
        
  TLorentzVector  lelectron( fElectron_MomX_Col, fElectron_MomY_Col, fElectron_MomZ_Col, fElectron_Energy_Col);

  return lelectron;	 

}

Double_t KPlus_Production::Get_Phi_X_LeptonPlane_RF () {

  fCos_Phi_X_LeptonPlane_RF = ( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ); // hep-ph/0410050v2
  fSin_Phi_X_LeptonPlane_RF = ( ( v3LxP.Dot( v3PhotonUnit  ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ); // hep-ph/0410050v2    
  if ( fSin_Phi_X_LeptonPlane_RF >= 0 )
    fPhi_X_LeptonPlane_RF   = fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) );
  if ( fSin_Phi_X_LeptonPlane_RF < 0 )
    fPhi_X_LeptonPlane_RF   = 360.0 - std::abs( fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ) );

  return fPhi_X_LeptonPlane_RF;

}

Double_t KPlus_Production::Get_Phi_TargPol_LeptonPlane_RF () {

  fCos_Phi_TargPol_LeptonPlane_RF = ( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ); // hep-ph/0410050v2
  fSin_Phi_TargPol_LeptonPlane_RF = ( ( v3LxS.Dot( v3PhotonUnit  ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ); // hep-ph/0410050v2
  if ( fSin_Phi_TargPol_LeptonPlane_RF >= 0 )
    fPhi_TargPol_LeptonPlane_RF = fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) );
  if ( fSin_Phi_TargPol_LeptonPlane_RF < 0 )
    fPhi_TargPol_LeptonPlane_RF = 360.0 - std::abs( fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ) );

  return fPhi_TargPol_LeptonPlane_RF;

}

Double_t KPlus_Production::Get_Total_Cross_Section() {

  Double_t total_sig;

  if(rParticle_scat_hadron == "Lambda"){
    total_sig = GetKPlus_CrossSection();
    //total_sig = GetKLambda_CrossSection(); -> SJDK 08/02/22 Not implemented yet
  }
  else if (rParticle_scat_hadron == "Sigma0"){
    total_sig = GetKPlus_CrossSection();
    //total_sig = GetKSigma_CrossSection(); -> SJDK 08/02/22 Not implemented yet
  }

  return total_sig;

}

/*--------------------------------------------------*/
/// Charged Pi+ module: 
/// Author: Z. Ahmed 
/// Date: 2017
/// Modified by A. Usman to switch to Kaon production - June/July 2021

Double_t  KPlus_Production::GetKPlus_CrossSection(){

  double_t sig_total;

  // --------------------------------------------------------------------------------------------------
  // CKY sigma L and T starts
  // --------------------------------------------------------------------------------------------------
  double lpar0 = 0., lpar1 = 0., lpar2 = 0., lpar3 = 0., lpar4 = 0., lpar5 = 0., lpar6 = 0.;
  double tpar0 = 0., tpar1 = 0., tpar2 = 0., tpar3 = 0., tpar4 = 0.;


  lpar0 = 0.;    lpar1 = 0.;    lpar2 = 0.;    lpar3 = 0.;    lpar4 = 0.;    lpar5 = 0.;    lpar6 = 0.;
  tpar0 = 0.;    tpar1 = 0.;    tpar2 = 0.;    tpar3 = 0.;    tpar4 = 0.;
 
  fSig_L = 0;
  fSig_T = 0;
 
  if ( ( fT_GeV > 0. ) && ( fT_GeV < 0.15 ) ) {
    eicSigmaL( fW_GeV,  fQsq_GeV, lpar0, lpar1, lpar2 , lpar3 , lpar4 , lpar5 , lpar6 );
    TF1 *fitCKYLonglandau = new TF1("sigmaL","landau", 0.0 , 0.15 );
    fitCKYLonglandau->FixParameter( 0 , lpar0 );
    fitCKYLonglandau->FixParameter( 1 , lpar1 );
    fitCKYLonglandau->FixParameter( 2 , lpar2 );
    fSig_L = fitCKYLonglandau->Eval(fT_GeV);
    if ( lpar0 == 0 || lpar1 == 0 || lpar2 == 0 )
      fSig_L = 0;
    fitCKYLonglandau = NULL;
    delete fitCKYLonglandau;
  }
  else if ( ( fT_GeV > 0.15 ) && ( fT_GeV < 0.5 ) ) {
    eicSigmaL( fW_GeV,  fQsq_GeV, lpar0, lpar1, lpar2 , lpar3 , lpar4 , lpar5 , lpar6 );
    TF1 *fitCKYLongexpo1 = new TF1("sigmaL","expo", 0.15 , 0.5 );
    fitCKYLongexpo1->FixParameter( 0 , lpar3 );
    fitCKYLongexpo1->FixParameter( 1 , lpar4 );
    fSig_L = fitCKYLongexpo1->Eval(fT_GeV);
    if ( lpar3 == 0 || lpar4 == 0 )
      fSig_L = 0;
    fitCKYLongexpo1 = NULL;
    delete fitCKYLongexpo1;
  }
  else if ( ( fT_GeV > 0.5 ) && ( fT_GeV < 1.3 ) ) {
    eicSigmaL( fW_GeV,  fQsq_GeV, lpar0, lpar1, lpar2 , lpar3 , lpar4 , lpar5 , lpar6 );
    TF1 *fitCKYLongexpo2 = new TF1("sigmaL","expo", 0.5 , 1.3 );
    fitCKYLongexpo2->FixParameter( 0 , lpar5 );
    fitCKYLongexpo2->FixParameter( 1 , lpar6 );
    fSig_L = fitCKYLongexpo2->Eval(fT_GeV);
    if ( lpar5 == 0 || lpar6 == 0 )
      fSig_L = 0;
    fitCKYLongexpo2 = NULL;
    delete fitCKYLongexpo2;
  }
  else {
    fSig_L = 0;
  }
   // SJDK - 02/06/22 - The validity range here was inconsistent, this only went from 0.0 to 0.15, leaving a gap between 0.15 to 0.2
  // I changed the range to remove this gap. 
  if ( ( fT_GeV > 0.0 ) && ( fT_GeV < 0.2 ) ) {
    eicSigmaT( fW_GeV,  fQsq_GeV, tpar0, tpar1, tpar2 , tpar3 , tpar4 );
    TF1 *fitCKYTranspol2 = new TF1("sigmaL","pol2", 0.0 , 0.2 );
    fitCKYTranspol2->FixParameter( 0 , tpar0 );
    fitCKYTranspol2->FixParameter( 1 , tpar1 );
    fitCKYTranspol2->FixParameter( 2 , tpar2 );
    fSig_T = fitCKYTranspol2->Eval(fT_GeV);
    if ( tpar0 == 0 || tpar1 == 0 || tpar2 == 0 )
      fSig_T = 0;
    fitCKYTranspol2 = NULL;
    delete fitCKYTranspol2;
  }
  else if ( ( fT_GeV > 0.2 ) && ( fT_GeV < 1.3 ) ) {
    eicSigmaT( fW_GeV,  fQsq_GeV, tpar0, tpar1, tpar2 , tpar3 , tpar4 );
    TF1 *fitCKYTransexpo = new TF1("sigmaL","expo", 0.2 , 1.3 );
    fitCKYTransexpo->FixParameter( 0 , tpar3 );
    fitCKYTransexpo->FixParameter( 1 , tpar4 );
    fSig_T = fitCKYTransexpo->Eval(fT_GeV);
    if ( tpar3 == 0 || tpar4 == 0 )
      fSig_T = 0;
    fitCKYTransexpo = NULL;
    delete fitCKYTransexpo;
  }
 
  // ------------------------------------------------------------------------------------------------
  // Improving the Kaon Sigma_L following GH's fortran code (pole_ration.f) located in main driectory
  // Added by AU on July 10, 2021
  // ------------------------------------------------------------------------------------------------

  double mkg = 0, mpig = 0, hbarc = 0, gpoleKL = 0, gpoleKS = 0, gpolepi = 0, gpoleKhyp = 0, r2_dip = 0, r2_mono = 0, Fpi_mono = 0, Fpi_dip = 0, pmono = 0, Fpi_fit = 0, fpisq = 0, Fkk = 0, fksq = 0, lNpi = 0, lNk = 0, gkhypn = 0, gpinn = 0, dl_poleKhyp = 0, dl_polepi = 0, ratio_Khyp_pi = 0, FT_GeV_neg = 0;

  gpoleKL = -13.3;
  gpoleKS = -3.5; 
  gpolepi = 13.1;

  r2_dip = 0.411;
  r2_mono = 0.431;                                                                                                                                                                                                              

  hbarc = 0.197;
  mkg = fX_Mass_GeV;
  mpig = 0.13957;

  FT_GeV_neg = -1.0 * fT_GeV;

  if ( recoil_hadron = Lambda ) {
    gpoleKhyp = gpoleKL;
  }
  else if ( recoil_hadron = Sigma0 ) {
    gpoleKhyp = gpoleKS;
  }

  Fpi_mono = 1.0 / ( 1.0 + ( r2_mono * fQsq_GeV) / (6 * (hbarc * hbarc)));   
  Fpi_dip = 1.0 / ( 1.0 + ( r2_dip * fQsq_GeV) / (12 * (hbarc * hbarc))) * ( 1.0 + ( r2_dip * fQsq_GeV) / (12 * (hbarc * hbarc)));
  pmono = 0.85;
  
  Fpi_fit = pmono * Fpi_mono + (1 - pmono) * Fpi_dip;
  fpisq = Fpi_fit * Fpi_fit;

  Fkk = 0.9 / (1.0 + fQsq_GeV / 0.462 );
  fksq = Fkk * Fkk;

  lNpi = 0.44;
  lNk = (0.44+0.80)/2;

  gkhypn = gpoleKhyp * ((lNk * lNk) - (mkg * mkg)) / ((lNk * lNk) - FT_GeV_neg); 
  gpinn = gpolepi * ((lNpi * lNpi) - (mpig * mpig)) / ((lNpi * lNpi) - FT_GeV_neg);

  dl_poleKhyp = ((gkhypn * gkhypn) * fksq) / ((FT_GeV_neg - (mkg * mkg))*(FT_GeV_neg - (mkg * mkg)));  
  dl_polepi = ((gpinn * gpinn) * fpisq) / ((FT_GeV_neg - (mpig * mpig))*(FT_GeV_neg - (mpig * mpig)));
  
  ratio_Khyp_pi = dl_poleKhyp / dl_polepi; 

  // --------------------------------------------------------------------------------

  fSig_VR = (0.1* fSig_T) + fEpsilon * (ratio_Khyp_pi* fSig_L);
  sig_total = fSig_VR;

  return sig_total;
}

/*--------------------------------------------------*/
/// Output generator detail

void KPlus_Production::Detail_Output() {

  ppiDetails << "Total events tried                                           " << setw(20) << fNGenerated   << endl;
  ppiDetails << "Total events recorded                                        " << setw(20) << fNRecorded    << endl;
  // 09/02/22 - Commented out, not used anymore
  //ppiDetails << "Max weight value                                             " << setw(20) << fEventWeightCeil << endl; 
  ppiDetails << "Number of events with w more than 10.6                       " << setw(20) << w_ev          << endl;
  ppiDetails << "Number of events with wsq negative                           " << setw(20) << w_neg_ev      << endl;
  ppiDetails << "Number of events with qsq less than 5                        " << setw(20) << qsq_ev        << endl;
  ppiDetails << "Number of events with Meson (X) energy NaN                   " << setw(20) << fNaN          << endl;
  ppiDetails << "Number of events failing conservation law check              " << setw(20) << fConserve     << endl;
  ppiDetails << "Total events passing conservation laws                       " << setw(20) << conserve   << endl;
  ppiDetails << "Total events failed energy conservation                      " << setw(20) << ene   << endl; 
  ppiDetails << "Total events failed momentum conserveation                   " << setw(20) << mom   << endl;

  ppiDetails << "Number of events with -t more than threshold                 " << setw(20) << t_ev          << endl;
  ppiDetails << "Number of events with w less than threshold                  " << setw(20) << fWSqNeg       << endl;
  ppiDetails << "Number of events with mom not conserve                       " << setw(20) << fNMomConserve << endl;
  ppiDetails << "Number of events with Sigma negative                         " << setw(20) << fNSigmaNeg    << endl;
  ppiDetails << "Number of lund events                                        " << setw(20) << fLundRecorded << endl;

  ppiDetails << "Seed used for the Random Number Generator                    " << setw(20) << fSeed         << endl;

}

/*--------------------------------------------------*/
/// Output format functions follow

void KPlus_Production::Lund_Output() {

  ppiOut << "3"
	 << " \t " << fPhi           // var 1
	 << " \t " << fPhiS          // var 2
	 << " \t " << fx             // var 3
	 << " \t " << "1"	       
	 << " \t " << fQsq_GeV       // var 4
	 << " \t " << fT_GeV         // var 5
	 << " \t " << fW_GeV 	       // var 6
	 << " \t " << fEpsilon       // var 7
	 << " \t " << fEventWeight   // var 8	   
	 << endl;
       
  // Produced Particle X
  ppiOut << setw(10) << "1" 
	 << setw(10) << "1" 
	 << setw(10) << "1" 
    // 	   << setw(10) << "11111111111"
	 << setw(10) << PDGtype(produced_X)
	 << setw(10) << "0" 
	 << setw(10) << "0" 
	 << setw(16) << r_lX_g.X()
	 << setw(16) << r_lX_g.Y()   
	 << setw(16) << r_lX_g.Z()  
	 << setw(16) << r_lX_g.E()
	 << setw(16) << fX_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;
     
  // Scattered electron
  ppiOut << setw(10) << "2" 
	 << setw(10) << "-1" 
	 << setw(10) << "1" 
	 << setw(10) << "11" 
	 << setw(10) << "0" 
	 << setw(10) << "0" 
	 << setw(16) << r_lscatelecg.X() 
	 << setw(16) << r_lscatelecg.Y() 
	 << setw(16) << r_lscatelecg.Z() 
	 << setw(16) << r_lscatelecg.E()
	 << setw(16) << fElectron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;
 	  
  // Recoiled hadron
  ppiOut << setw(10) << "3" 
	 << setw(10) << "1" 
	 << setw(10) << "1" 
	 << setw(10) << PDGtype(recoil_hadron)
	 << setw(10) << "0" 
	 << setw(10) << "0" 
	 << setw(16) << r_l_scat_hadron_g.X() 
	 << setw(16) << r_l_scat_hadron_g.Y()
	 << setw(16) << r_l_scat_hadron_g.Z()
	 << setw(16) << r_l_scat_hadron_g.E()
	 << setw(16) << f_Scat_hadron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;
}

void KPlus_Production::KPlus_Pythia6_Out_Init() {

	print_itt = 0;

	ppiOut << "DEMP Event FILE" << endl;
	ppiOut << "============================================" << endl;
	ppiOut << "I, ievent, nParticles, Weight" << endl;
	ppiOut << "============================================" << endl;
	ppiOut << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)" << endl;
	ppiOut << "============================================" << endl;

}

void KPlus_Production::KPlus_Pythia6_Output() {

  ppiOut << "0" << " \t\t\t "  << print_itt << " \t\t\t " << "1" << " \t\t\t " << fEventWeight << endl;           // var 1

	print_itt++;

	ppiOut << "============================================" << endl;

 	///*--------------------------------------------------*/
  	// Initial State
 
      ppiOut  << "1" 
  	   << setw(6) << "21" 
  	   << setw(6) << "11"
  	   << setw(6) << "0" 
  	   << setw(6) << "3" 
  	   << setw(6) << "4" 

  	   << setw(14) << r_lelectrong.X()
  	   << setw(14) << r_lelectrong.Y()   
  	   << setw(14) << r_lelectrong.Z()
  	   << setw(14) << r_lelectrong.E()
  	   << setw(14) << fElectron_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

      ppiOut << "2" 
  	   << setw(6) << "21" 
  	   << setw(6) << "2212"
  	   << setw(6) << "0" 
  	   << setw(6) << "5" 
  	   << setw(6) << "6" 

  	   << setw(14) << r_lprotong.X()
  	   << setw(14) << r_lprotong.Y()   
  	   << setw(14) << r_lprotong.Z()  
  	   << setw(14) << r_lprotong.E()
  	   << setw(14) << fProton_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

      ppiOut << "3" 
  	   << setw(6) << "21" 
  	   << setw(6) << "22"
  	   << setw(6) << "1" 
  	   << setw(6) << "0" 
  	   << setw(6) << "0" 

  	   << setw(14) << r_lphotong.X()
  	   << setw(14) << r_lphotong.Y()   
  	   << setw(14) << r_lphotong.Z()  
  	   << setw(14) << r_lphotong.E()
  	   << setw(14) << r_lphotong.M()
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

 	///*--------------------------------------------------*/
  	// Final State
      
      // Scattered electron
      ppiOut << "4" 
  	   << setw(6) << "1" 
  	   << setw(6) << "11" 
  	   << setw(6) << "1" 
  	   << setw(6) << "0"
  	   << setw(6) << "0"
 
  	   << setw(14) << r_lscatelecg.X() 
  	   << setw(14) << r_lscatelecg.Y() 
  	   << setw(14) << r_lscatelecg.Z() 
  	   << setw(14) << r_lscatelecg.E()
  	   << setw(14) << fElectron_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;
  	  
      // Recoiled hadron
      ppiOut << "5" 
  	   << setw(6) << "1" 
  	   << setw(6) << PDGtype(recoil_hadron)
  	   << setw(6) << "2" 
  	   << setw(6) << "0"
  	   << setw(6) << "0"
 
  	   << setw(14) << r_l_scat_hadron_g.X() 
  	   << setw(14) << r_l_scat_hadron_g.Y()
  	   << setw(14) << r_l_scat_hadron_g.Z()
  	   << setw(14) << r_l_scat_hadron_g.E()
  	   << setw(14) << f_Scat_hadron_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;
 
      // Produced Particle X
      ppiOut << "6" 
  	   << setw(6) << "1" 
  	   << setw(6) << PDGtype(produced_X)
  	   << setw(6) << "2" 
  	   << setw(6) << "0" 
  	   << setw(6) << "0"

  	   << setw(14) << r_lX_g.X()
  	   << setw(14) << r_lX_g.Y()   
  	   << setw(14) << r_lX_g.Z()  
  	   << setw(14) << r_lX_g.E()
  	   << setw(14) << fX_Mass_GeV
  	   << setw(6) << fVertex_X
  	   << setw(6) << fVertex_Y
  	   << setw(6) << fVertex_Z
  	   << endl;

	ppiOut << "=============== Event finished ===============" << endl;

}

void KPlus_Production::KPlus_HEPMC3_Out_Init() {
 
  print_itt = 0;
  ppiOut << "HepMC::Version 3.02.02" << endl;
  ppiOut << "HepMC::Asciiv3-START_EVENT_LISTING" << endl;

}

void KPlus_Production::KPlus_HEPMC3_Output() {
  
  // HEPMC3 output for Athena simulations

  // First line - E - Event# - #Vertices - #Particles
  ppiOut << "E" << " "  << print_itt <<  " " << "1" << " " << 5 << endl;
  print_itt++;
  // Second line, Units - U - ENERGY UNIT - DISTANCE UNIT
  ppiOut << "U" << " " << "GEV" << " " << "MM" << endl;
  // Third line, optional attributes, the weight
  ppiOut << "A" << " " << "0" << " " << "weight" << " " <<  fEventWeight << endl;
  // Beam particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (4, incoming beam particle)
  ppiOut << "P" << " " << "1" << " " << "0" << " " << "11" << " " << r_lelectrong.X() << " " << r_lelectrong.Y() << " " << r_lelectrong.Z() << " " << r_lelectrong.E() << " " << fElectron_Mass_GeV << " " << "4" << endl;
  ppiOut << "P" << " " << "2" << " " << "0" << " " << "2212" << " " << r_lprotong.X() << " " << r_lprotong.Y() << " " << r_lprotong.Z() << " " << r_lprotong.E() << " " << fProton_Mass_GeV << " " << "4" << endl;
  // Vertex line - V - 1 - 0 - [1,2]
  ppiOut << "V" << " " << "-1" << " " << "0" << " " << "[1,2]" << endl;
  // Output particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (1, undecayed physical particle)
  // Scattered electron
  ppiOut << "P" << " " << "3" << " " << "-1" << " " << "11" << " " << r_lscatelecg.X() << " "  << r_lscatelecg.Y() << " "  << r_lscatelecg.Z() << " " << r_lscatelecg.E() << " " << fElectron_Mass_GeV << " " << "1" << endl;
  // Produced meson
  ppiOut << "P" << " " << "4" << " " << "-1" << " " << PDGtype(produced_X) << " " << r_lX_g.X() << " "  << r_lX_g.Y() << " "  << r_lX_g.Z() << " " << r_lX_g.E() << " " << fX_Mass_GeV << " " << "1" << endl;
  // Recoil nucleon
  ppiOut << "P" << " " << "5" << " " << "-1" << " " << PDGtype(recoil_hadron) << " " << r_l_scat_hadron_g.X() << " "  << r_l_scat_hadron_g.Y() << " "  << r_l_scat_hadron_g.Z() << " " << r_l_scat_hadron_g.E() << " " << f_Scat_hadron_Mass_GeV << " " << "1" << endl;
  
}
