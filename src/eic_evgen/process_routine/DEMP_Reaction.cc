#include "reaction_routine.h"
#include "eic.h"
// Love Preet 27/01/23 - Added TF1 include here for now as it is needed in the KPlus function
#include "TF1.h"

using namespace std;

DEMP_Reaction::DEMP_Reaction() { 

  cout << "Program Start" << endl;

}

/*--------------------------------------------------*/
/// DEMP_Reaction 

DEMP_Reaction::DEMP_Reaction(TString particle_str, TString hadron_str) { 

  rParticle = particle_str;
  rHadron = hadron_str;

}

DEMP_Reaction::~DEMP_Reaction() {

  DEMPOut.close();
  DEMPDetails.close();

}

void DEMP_Reaction::process_reaction() {
 
  Init();
  
  if (gOutputType == "Pythia6"){
    DEMPReact_Pythia6_Out_Init();
  }
  else if (gOutputType == "HEPMC3"){
    DEMPReact_HEPMC3_Out_Init();
  }

  for( long long int i = 0; i < rNEvents; i++ ) {
 
    rNEvent_itt = i;
    fNGenerated ++;
 
    Progress_Report();  // This is happens at each 10% of the total event is processed
    Processing_Event();
  }
 
  Detail_Output();
 
}

void DEMP_Reaction::Init() {

  pim* myPim;

  pd = dynamic_cast<pim*>(myPim);
	
  rParticle_charge = ExtractCharge(rParticle);

  sTFile = Form("./LundFiles/eic_%s.txt", gfile_name.Data());
  sLFile= Form("./LundFiles/eic_input_%s.dat", gfile_name.Data());
   
  DEMPOut.open( sLFile.c_str() );
  DEMPDetails.open( sTFile.c_str() );
	
  qsq_ev = 0, t_ev = 0, w_neg_ev = 0, w_ev = 0;
  rNEvents = fNEvents;
  rNEvent_itt = 0;

  // 02/06/21 SJDK 
  // Set these values once the beam energies are read in
  fPSF = ( fEBeam * ( fScatElec_E_Hi - fScatElec_E_Lo ) *( sin( fScatElec_Theta_F ) - sin( fScatElec_Theta_I ) ) * 2 * fPI *( sin( fEjectileX_Theta_F ) - sin( fEjectileX_Theta_I ) ) * 2 * fPI );
  fElectron_Kin_Col_GeV = fEBeam;
  fElectron_Kin_Col = fElectron_Kin_Col_GeV * 1000.0;

  // cout << rNEvents << "    " << fNEvents << endl;
	
  rFermiMomentum = pd->fermiMomentum();

  // ----------------------------------------------------
  // Proton in collider (lab) frame

  r_lproton = GetProtonVector_lab();
  r_lprotong = GetProtonVector_lab() * fm;

  // ----------------------------------------------------
  // Electron in collider (lab) frame

  cout << "Fermi momentum: " << rFermiMomentum << endl;

  r_lelectron	 = GetElectronVector_lab();
  r_lelectrong = r_lelectron * fm;

  ///*--------------------------------------------------*/
  /// Getting the particle mass from the data base
 
  produced_X = ParticleEnum(rParticle);
  fX_Mass = ParticleMass(produced_X)*1000; //MeV
  fX_Mass_GeV = fX_Mass/1000; //GeV

  cout << rParticle << "  " << produced_X << "  " << fX_Mass_GeV <<  endl;
  cout << rParticle_charge << endl;

  if (rHadron == "Neutron" ) {
    rParticle_scat_hadron  = "Neutron"; 
    recoil_hadron  = Neutron; 
    f_Scat_hadron_Mass     = fNeutron_Mass;
    f_Scat_hadron_Mass_GeV = f_Scat_hadron_Mass/1000;
  }
  else if (rHadron == "Proton" ) {	
    rParticle_scat_hadron  = "Proton"; 
    recoil_hadron  = Proton;
    f_Scat_hadron_Mass = fProton_Mass;
    f_Scat_hadron_Mass_GeV = f_Scat_hadron_Mass/1000;
  } 
  else if(rHadron == "Lambda"){
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

  rDEG2RAD   = fPI/180.0;

  fX_Theta_I = fEjectileX_Theta_I ;
  fX_Theta_F = fEjectileX_Theta_F;

  cout << "Produced particle in exclusive production: " << rParticle << ";  with mass: " << fX_Mass << " MeV "<< endl;
  cout << fEBeam << " GeV electrons on " << fPBeam << " GeV ions" << endl;
  
  // Set luminosity value based upon beam energy combination, note that if no case matches, a default of 1e33 is assumed. Cases are a set of nominal planned beam energy combinations for the EIC (and EICC)
  // See slide 11 in https://indico.cern.ch/event/1072579/contributions/4796856/attachments/2456676/4210776/CAP-EIC-June-7-2022-Seryi-r2.pdf
  // If available in the future, this could be replaced by some fixed function
  if ((fEBeam == 5.0 ) && (fPBeam == 41.0) ){
    fLumi = 0.44e33;
  }
  else if ((fEBeam == 5.0 ) && (fPBeam == 100.0) ){
    fLumi = 3.68e33;
  }
  else if ((fEBeam == 10.0 ) && (fPBeam == 100.0) ){
    fLumi = 4.48e33;
  }
  else if ((fEBeam == 18.0 ) && (fPBeam == 275.0) ){
    fLumi = 1.54e33;
  }
  else if ((fEBeam == 3.5 ) && (fPBeam == 20) ){ // EICC optimal beam energy combination
    fLumi = 2e33;
  }
  else if ((fEBeam == 2.8 ) && (fPBeam == 13) ){ // EICC lowest beam energy combination
    fLumi = 0.7e33;
  }
  else{
    cout << "!!! Notice !!! The beam energy combination simulated does not match an expected case, a default luminosity value of - " << fLumi << " cm^2s^-1 has been assumed. !!! Notice !!!" << endl;
  }
  
}

void DEMP_Reaction::Processing_Event() {

  // ----------------------------------------------------
  // Considering Fermi momentum for the proton
  // ----------------------------------------------------
  // SJDK - 31/01/23 - This doesn't seem to do anything?
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
  // SJDK 31/01/23 - Again, should this be a user defined/set range?
  // SJDK 30/01/23 - Changed Qsq range to match new validity range from Love's paramaterisation
  if ( fQsq_GeV < 1.0 || fQsq_GeV > 35.0 ) {
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

  // SJDK 31/01/23 - Should have the lower/upper W ranges be user defined?
  // Also, moved this cut up in the code, W is known at this point, so why wait to cut on it?
  //if ( fW_GeV < 3.0 || fW_GeV > 10.6 ) { // SJDK 31/01/23 - Previous range utilised
  if ( fW_GeV < 2.0 || fW_GeV > 10 ) { // SJDK 31/01/23 - New range for W to work with K+ cross section model
    w_ev++;
    return;
  }
  
  // 13/12/22 - SJDK - This is the start of the block that will need to be replaced by the ROOT function Rory used to determine the pion momentum
  // 21/12/22 - SJDK - Should split this out into its own class, then have two different variants (Rory vs Ahmed)
  // ---------------------------------------------------------
  // Pion momentum in collider frame, analytic solution starts
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
  // SJDK 27/01/23 - For Kaon events, 0.5 is too stringent for conservation law check with the current (Ahmed) method to determine meson properties
  // Hopefully, Rory's method will be better here and we can utilise the same conservation law check for both particles
  if (rParticle == "Pi+"){
    if ( pd->CheckLaws( r_lelectron, r_lproton, r_lscatelec, r_lX, r_l_scat_hadron, 0.5) != 1 ){
      fConserve++;
      return;
    }
    else if (rParticle == "K+"){
      if ( pd->CheckLaws( r_lelectron, r_lproton, r_lscatelec, r_lX, r_l_scat_hadron, 10) != 1 ){
	fConserve++;
	return;
      }
    }
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
  
  // 31/01/23 - SJDK - Kinematics type 1 is FF and 2 is TSSA, this reaction class shouldn't care about this and only have a limit depending upon particle type?
  /*
  if ( gKinematics_type == 1 && fT_GeV > 0.5 ) {
    t_ev++;
    return;
  }
     
  if ( gKinematics_type == 2 && fT_GeV > 1.3 ) {
    t_ev++;
    return;
  }
  */ 

  // 31/01/23 SJDK - New limit on t, remove only events outside the parameterisation range, limits depend upon particle type, need to add Pi0 version
  if (rParticle == "Pi+" && fT_GeV > 1.3 ) {
    t_ev++;
    return;
  }
  else if (rParticle == "K+" && fT_GeV > 2.0) {
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

  /*--------------------------------------------------*/
  // Get the Phi scattering angle with respect to the electron scattering plane
  fPhi   = Get_Phi_X_LeptonPlane_RF ();

  /*--------------------------------------------------*/
  // Get the Phi scattering angle with respect to the electron scattering plane
  fPhiS  = Get_Phi_TargPol_LeptonPlane_RF();

  fTheta_X_Photon_RF       = fRAD2DEG * acos( ( v3Photon.Dot( v3X     ) ) / ( v3Photon.Mag()  * v3X.Mag()    ) );
  if ( fTheta_X_Photon_RF < 0 ) { fTheta_X_Photon_RF = 180.0 + fTheta_X_Photon_RF; }

  // -----------------------------------------------------------------------------------
  // If we have fermi momentum then epsilon should be in rest frame 
  // The theta angle of scattered angle used in expression of epsilon is the angle 
  // with respect to direction of incoming electron in the rest frame of target hadron
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
  fJacobian_CM = ( (lphoton_rfg.Vect()).Mag() - fBeta_CM_RF * lphoton_rfg.E() ) / ( fGamma_CM_RF * ( 1.0 - pow(fBeta_CM_RF,2) ) ); // Eqn 22 in paper
 
  fA = fJacobian_CM * fX_Mom_CM_GeV / fPi; // Eqn 21 in paper
 
  // ----------------------------------------------------
  // Jacobian dOmega* / dOmega dimensionless
  // ----------------------------------------------------
  fJacobian_CM_RF  = ( pow((lX_rf.Vect()).Mag(),2)*fW) / 
    ( fX_Mom_CM * std::abs( ( fProton_Mass + lphoton_rf.E()) * (lX_rf.Vect()).Mag() - 
			    ( lX_rf.E() * (lphoton_rf.Vect()).Mag() * cos( lX_rf.Theta() ) ) ) ); // Differs from next line in photon vect -> lphoton_rf vs r_lphoton
 
  fJacobian_CM_Col = ( ( pow((r_lX.Vect()).Mag(),2) * fW ) / // This one is actually used subsequently, so this must be Eqn 20
		       ( fX_Mom_CM * std::abs( ( fProton_Mass + r_lphoton.E() ) * (r_lX.Vect()).Mag() -
					       ( r_lX.E() * (r_lphoton.Vect()).Mag() * cos( r_lX.Theta() ) ) ) ) ); 


  //	 cout <<  lX_rf.Vect().Mag() << "  " << << << << << << << << endl;
  //	 cout << fJacobian_CM_RF << "    " << fJacobian_CM_Col << endl;

  // -----------------------------------------------------------------------------------------------------------
  // CKY sigma L and T starts
  // -----------------------------------------------------------------------------------------------------------
  //	 r_fSig_T = 1;
  //	 r_fSig_L = 1;
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
  //             Lab cross section     Phase Space   Conversion     Luminosity                Total events tried
  // Hz        = ub / ( sr^2 * GeV ) * GeV * sr^2 * ( cm^2 / ub ) * ( # / ( cm^2 * sec ) ) / ( # )

  // SJDK 24/06/21 - Explicitly taking the absolute value of the weight such that the value is positive! Shouldn't matter since any -ve cross section events should be dumped above
  fEventWeight = abs(fSigma_Col * fPSF * fuBcm2 * fLumi / fNEvents);   // in Hz

  fNRecorded++;
  fLundRecorded++;
  fRatio = fNRecorded / fNGenerated;

  if (gOutputType == "Pythia6"){
      DEMPReact_Pythia6_Output();
  }
  else if (gOutputType == "LUND"){
    Lund_Output();
  }
  else if (gOutputType == "HEPMC3"){
    DEMPReact_HEPMC3_Output();
  }
}

void DEMP_Reaction::Progress_Report() {

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

TLorentzVector DEMP_Reaction::GetProtonVector_lab() {

  // Crossing angle
  //	 fProton_Theta_Col = 0.050;
  //	 fProton_Theta_Col = 0.025;
  // SJDK - 12/01/22
  // Set crossing angle to 0 for fun4all, also required for ATHENA simulations
  fProton_Theta_Col = 0.0;

  ///*--------------------------------------------------*/
  /// The 
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

void DEMP_Reaction::Consider_Proton_Fermi_Momentum() {

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

TLorentzVector DEMP_Reaction::GetElectronVector_lab() {

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

Double_t DEMP_Reaction::Get_Phi_X_LeptonPlane_RF () {

  fCos_Phi_X_LeptonPlane_RF = ( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ); // hep-ph/0410050v2
  fSin_Phi_X_LeptonPlane_RF = ( ( v3LxP.Dot( v3PhotonUnit  ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ); // hep-ph/0410050v2    
  if ( fSin_Phi_X_LeptonPlane_RF >= 0 )
    fPhi_X_LeptonPlane_RF   = fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) );
  if ( fSin_Phi_X_LeptonPlane_RF < 0 )
    fPhi_X_LeptonPlane_RF   = 360.0 - std::abs( fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ) );

  return fPhi_X_LeptonPlane_RF;

}

Double_t DEMP_Reaction::Get_Phi_TargPol_LeptonPlane_RF () {

  fCos_Phi_TargPol_LeptonPlane_RF = ( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ); // hep-ph/0410050v2
  fSin_Phi_TargPol_LeptonPlane_RF = ( ( v3LxS.Dot( v3PhotonUnit  ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ); // hep-ph/0410050v2
  if ( fSin_Phi_TargPol_LeptonPlane_RF >= 0 )
    fPhi_TargPol_LeptonPlane_RF = fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) );
  if ( fSin_Phi_TargPol_LeptonPlane_RF < 0 )
    fPhi_TargPol_LeptonPlane_RF = 360.0 - std::abs( fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ) );

  return fPhi_TargPol_LeptonPlane_RF;

}

Double_t DEMP_Reaction::Get_Total_Cross_Section() {

  Double_t total_sig, total_sig2;

  if (rParticle == "Pi+"){
    total_sig = GetPiPlus_CrossSection(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon);
  }
  else if (rParticle == "Pi0"){
    total_sig = GetPi0_CrossSection();
  }
  else if (rParticle == "K+"){
    total_sig = GetKPlus_CrossSection(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon, rHadron);
    // SJDK - 31/01/23 - This second case gets the total cross section for the K+ as calculated from Ali's earlier scaling calculation, retain this for comparison.
    total_sig2 = GetKPlus_CrossSection_Scaling(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon, fX_Mass_GeV, rHadron);
  }
  
  //cout << fT_GeV <<  "  " << fW_GeV << "   " << fQsq_GeV << "   " << "KPlus Paramterisation - " << total_sig << " !!! Old Scaling method - " << total_sig2 << endl;
  
  return total_sig;

}

// Love Preet 27/01/23 - Will implement this here for now and move it to a separate file later once it is tested and working
Double_t DEMP_Reaction::GetKPlus_CrossSection(double ft, double fw, double fqsq, double feps, TString fHadron) {

  double_t sig_total;
     
  double w,q2,t; // For the sake of testing, I redescribed the fw,fqsq,ft as w,q2,t as I alraedy had a file in my PC to calculate the crossection with these notations
  w=fw; q2=fqsq; t=ft;
     
  int w_1,w_2,q2_1,q2_2;// w_1, rounded value of the input W, w_2 is the next (or previous) value of W. q2_1 is the rounded input value of Q2, q2_2 is the next (or previous) value of Q2 in the array
  w_1 = round(w);  
  q2_1 = round(q2);
  //....................................................................................................................................................................
  // Below are set of if and else if conditions for boundary values of w and q2 
  //....................................................................................................................................................................
    
  if (w == w_1){
    w_1 = w;
    w_2 = w;
    
    if (q2 > q2_1){
      q2_2 = q2_1 + 1;
    }
    
    else if (q2 < q2_1){
      q2_2 = q2_1;
      q2_1 = q2_2 -1;
    }
     
    else if (q2 == q2_1 ){
      q2_1 = q2;
      q2_2 = q2;}
     
  }
   
  //................................................................................. 
   
  else if (w > w_1){
    w_2 = w_1 + 1;
    
    if (q2 > q2_1){
      q2_2 = q2_1 + 1;
    }
    
    else if (q2 < q2_1){
      q2_2 = q2_1;
      q2_1 = q2_2 -1;
    }
     
    else if (q2 == q2_1 ){
      q2_1 = q2;
      q2_2 = q2;}
  }
    
  //................................................................................... 
  else if (w < w_1){
    w_2 =w_1;
    w_1 = w_2 - 1;
      
    if (q2 > q2_1){
      q2_2 = q2_1 + 1;
    }
    
    else if (q2 < q2_1){
      q2_2 = q2_1;
      q2_1 = q2_2 -1;
    }
     
    else if (q2 == q2_1 ){
      q2_1 = q2;
      q2_2 = q2;}
  }

  //.................................................................................................................................................................... 
  // Calcualtion of t-min values at each w and q2 (i.e. at all the four combinations of w's and q2's ) 
  //.................................................................................................................................................................... 
  // t1 is the t_min value UNLESS the first fit fails
  // If the first fit fails, t2 is the t_min value. If the first fit did NOT fail, t2 is the intersection point between fits 1 and 2
  // If the second fit fails, t3 is the t_min value. If the second fit did NOT fail, t3 is the intersection point between fits 2 and 3
  
  int x_1,x_2,y_1,y_2; // Redefination of w's and q2's as the array start from 0 but w starts from 2 and the q2 starts from 1 in the file.
  x_1 = w_1-2;
  x_2 = w_2-2;
  y_1 = q2_1-1;
  y_2 = q2_2-1;  // Till this point variables are same for sigT and sigL

  //....................................................................................................................................................................     
  // Variables that need to be defined globally somewhere ( These variables are only for sigT and I have used different variables for sigL)
  double t_1,t_2,t_3,t_4; // 4 corners (points) in t to interpolate over
  double sigT1,sigT2,sigT3,sigT4; // Value of SigT at each corner of the square to interpolate over
  double lsigT1,lsigT2,lsigT3,lsigT4, slb,slt,sll,slr,altb,allr,fsigTLa,fsigTa; // Logarithm of sigma values and slopes along the bottom, top, left, right and the average of the slope along the top/bottom and left/right. fsigTLa is the result in logarithmic form and fsigTa is the result
  double sigT11,sigT14; // Bottom left and top right sigma points, these are calculated if there are less than 3 (or 2) points to interpolate over
  //...................................................................................................................................................................
  //................................................................................... 
 
  /*if(SigPar[1][x_1][y_1][10] != -10001){ // t_1 is t_min for the first corner - this is always the bottom left corner
    t_1 = SigPar[1][x_1][y_1][10];
    }
    else if (SigPar[1][x_1][y_1][11] != -10001){
    t_1 =  SigPar[1][x_1][y_1][11];
    } 
    else if ( SigPar[1][x_1][y_1][12] != -10001){
    t_1 =  SigPar[1][x_1][y_1][12];
    }
    //................................................................................... 
  
    if(SigPar[1][x_1][y_2][10] != -10001) {  // t_2 is t_min for the second corner - this is always the bottom right corner
    t_2 = SigPar[1][x_1][y_2][10];
    }
    else if (SigPar[1][x_1][y_2][11] != -10001){
    t_2 = SigPar[1][x_1][y_2][11];
    }
    else if (SigPar[1][x_1][y_2][12] != -10001){
    t_2 = SigPar[1][x_1][y_2][12];
    }*/
  //...................................................................................        
 
  if(SigPar[1][x_2][y_1][10] != -10001){ // t_3 is t_min for the third corner - this is always the top left corner
    t_3 = SigPar[1][x_2][y_1][10];
  }
  else if (SigPar[1][x_2][y_1][11] != -10001){
    t_3 = SigPar[1][x_2][y_1][11];
  }
  else if (SigPar[1][x_2][y_1][12] != -10001){
    t_3 = SigPar[1][x_2][y_1][12];
  }
  else {
    return -100;
  }
  //...................................................................................
  
  /* if(SigPar[1][x_2][y_2][10] != -10001){  // t_4 is t_min for the fourth corner - this is always the top right corner
     t_4 = SigPar[1][x_2][y_2][10];
     }
     else if (SigPar[1][x_2][y_2][11] != -10001){
     t_4 =SigPar[1][x_2][y_2][11];
     }
     else if (SigPar[1][x_2][y_2][12] != -10001){
     t_4 = SigPar[1][x_2][y_2][12];
     } */
  //.................................................................................................................................................................... 
  // Calcualtion of sigT's at all the corners of the square 
  //....................................................................................................................................................................  
 
  if (t>= t_3 && t<2.0){ // t_3 corner is the one that we will loose at the end

    //Calculating the sigT1 at bottom left corner of the square..................................................................................
    if (t>=SigPar[1][x_1][y_1][10] && t<SigPar[1][x_1][y_1][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_1][2]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_1][3]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_1][4]);
      sigT1=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_1][y_1][11] && t<SigPar[1][x_1][y_1][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_1][5]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_1][6]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_1][7]);
      sigT1=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_1][y_1][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_1][y_1][8]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_1][9]);
      sigT1=parasigT->Eval(t); 
	 
      if (SigPar[1][x_1][y_1][8] == 0 && SigPar[1][x_1][y_1][9] == 0){
	sigT1=0;
      }
    }
	 
    else {
      sigT1=0;
    }

    //Calculating the sigT2 at bottom right corner of the square..................................................................................
    if (t>=SigPar[1][x_1][y_2][10] && t<SigPar[1][x_1][y_2][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_2][2]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_2][3]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_2][4]);
      sigT2=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_1][y_2][11] && t<SigPar[1][x_1][y_2][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_2][5]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_2][6]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_2][7]);
      sigT2=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_1][y_2][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_1][y_2][8]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_2][9]);
      sigT2=parasigT->Eval(t); 
	 
      if (SigPar[1][x_1][y_2][8] == 0 && SigPar[1][x_1][y_2][9] == 0){
	sigT2=0;
      }
    }
	 
    else {
      sigT2=0;
    }

    //Calculating the sigT3 at top left corner of the square..................................................................................
    if (t>=SigPar[1][x_2][y_1][10] && t<SigPar[1][x_2][y_1][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_1][2]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_1][3]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_1][4]);
      sigT3=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_2][y_1][11] && t<SigPar[1][x_2][y_1][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_1][5]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_1][6]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_1][7]);
      sigT3=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_2][y_1][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_2][y_1][8]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_1][9]);
      sigT3=parasigT->Eval(t); 
	 
      if (SigPar[1][x_2][y_1][8] == 0 && SigPar[1][x_2][y_1][9] == 0){
	sigT3=0;
      }
    }
	 
    else {
      sigT3=0;
    }

    //Calculating the sigT4 at top right corner of the square..................................................................................

    if (t>=SigPar[1][x_2][y_2][10] && t<SigPar[1][x_2][y_2][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_2][2]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_2][3]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_2][4]);
      sigT4=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_2][y_2][11] && t<SigPar[1][x_2][y_2][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_2][5]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_2][6]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_2][7]);
      sigT4=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_2][y_2][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_2][y_2][8]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_2][9]);
      sigT4=parasigT->Eval(t); 
	 
      if (SigPar[1][x_2][y_2][8] == 0 && SigPar[1][x_2][y_2][9] == 0){
	sigT4=0;
      }
    }
	 
    else {
      sigT4=0;
    }
    //....................................................................................................................................................................
    // Difterent if and else conditions to find the crossection values at the given points.
    //....................................................................................................................................................................
    // SJDK - 31/01/23 - This seems to spam things to screen, even when running pi+
    /*
      if (sigT1 == sigT2 && sigT2 == sigT3 && sigT3 == sigT4 && sigT4 == sigT1){ // if the w and q2 will have whole number values
      cerr<<"fsigTa = "<<sigT1 <<endl;
      }   
    */
    //...................................................................................

    if (sigT1 != 0 && sigT2 != 0 && sigT3 != 0 && sigT4 != 0){ // if all the four corners are  present
   
      // Taking the log of claculated sigT values    
      lsigT1 = TMath::Log(sigT1); //log value of sigT1.
      lsigT2 = TMath::Log(sigT2); //log value of sigT2.
      lsigT3 = TMath::Log(sigT3); //log value of sigT3.
      lsigT4 = TMath::Log(sigT4); //log value of sigT4.
      // Calculated slopes of different lines
      slb = lsigT4-lsigT3; //interpolation from the third corner (i.e. top left corner)
      slt = lsigT2-lsigT1;
      sll = -(lsigT1-lsigT3);
      slr = -(lsigT2-lsigT4);
      // Taking averages of the slopes
      altb = ((slb +slt)/2 );
      allr = ((sll +slr)/2 );
        
      //  Applying taylor's series formula for the average slopes
      fsigTLa = lsigT3 + (q2-q2_1)*altb + (w-w_2)*allr;
      
      // Find the anti-log of the taylor's series formula value  
      fsigTa = exp(fsigTLa);
    }
    //...................................................................................

    else if (sigT1 == 0 && sigT2 == 0 && sigT4 == 0)  // if we loose the 1st, 2nd and the 4th corner simultaneously 
      {
	// In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) will always be there and for other two corners, find the value of the cross-section at the first and the fourth corner at the minimum value of t. After that we can interpolate them.

	// First try to find t_1 and t_4 

	if(SigPar[1][x_1][y_1][10] != -10001){ // t_1 is t_min for the first corner - this is always the bottom left corner
	  t_1 = SigPar[1][x_1][y_1][10];
	}
	else if (SigPar[1][x_1][y_1][11] != -10001){
	  t_1 =  SigPar[1][x_1][y_1][11];
	} 
	else if ( SigPar[1][x_1][y_1][12] != -10001){
	  t_1 =  SigPar[1][x_1][y_1][12];
	}
	else {
	  return -100;
	}
 
	//...................................................................................  
	if(SigPar[1][x_2][y_2][10] != -10001){  // t_4 is t_min for the fourth corner - this is always the top right corner
	  t_4 = SigPar[1][x_2][y_2][10];
	}
	else if (SigPar[1][x_2][y_2][11] != -10001){
	  t_4 =SigPar[1][x_2][y_2][11];
	}
	else if (SigPar[1][x_2][y_2][12] != -10001){
	  t_4 = SigPar[1][x_2][y_2][12];
	}
	else {
	  return -100;
	}


	//Calculating the sigT11 at bottom left corner of the square      
	if (t_1>=SigPar[1][x_1][y_1][10]&& t_1<SigPar[1][x_1][y_1][11]){ 
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0, SigPar[1][x_1][y_1][2]); 
	  parasigT->FixParameter(1, SigPar[1][x_1][y_1][3]);
	  parasigT->FixParameter(2, SigPar[1][x_1][y_1][4]);
	  sigT11=parasigT->Eval(t_1);
	}
	 
	else if (t_1>=SigPar[1][x_1][y_1][11] && t_1<SigPar[1][x_1][y_1][12]){
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0,SigPar[1][x_1][y_1][5]); 
	  parasigT->FixParameter(1,SigPar[1][x_1][y_1][6]);
	  parasigT->FixParameter(2,SigPar[1][x_1][y_1][7]);
	  sigT11=parasigT->Eval(t_1); 
	} 
	
	else if (t_1>=SigPar[1][x_1][y_1][12] && t_1<2.0){
	  TF1* parasigT= new TF1("parasigT","expo");
	  parasigT->FixParameter(0, SigPar[1][x_1][y_1][8]); 
	  parasigT->FixParameter(1, SigPar[1][x_1][y_1][9]);
	  sigT11=parasigT->Eval(t_1); 
	  
	  if (SigPar[1][x_1][y_1][8] == 0 && SigPar[1][x_1][y_1][9] == 0){
	    sigT11 =0;
	  }
	}
	 
	else {
	  sigT11=0;
	}

	//Calculating the sigT14 at bottom left corner of the square
   
	if (t_4>=SigPar[1][x_2][y_2][10]&& t_4<SigPar[1][x_2][y_2][11]){ 
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0, SigPar[1][x_2][y_2][2]); 
	  parasigT->FixParameter(1, SigPar[1][x_2][y_2][3]);
	  parasigT->FixParameter(2, SigPar[1][x_2][y_2][4]);
	  sigT14=parasigT->Eval(t_4);
	}
	 
	else if (t_4>=SigPar[1][x_2][y_2][11] && t_1<SigPar[1][x_2][y_2][12]){
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0,SigPar[1][x_2][y_2][5]); 
	  parasigT->FixParameter(1,SigPar[1][x_2][y_2][6]);
	  parasigT->FixParameter(2,SigPar[1][x_2][y_2][7]);
	  sigT14=parasigT->Eval(t_4); 
	} 
	
	else if (t_4>=SigPar[1][x_1][y_1][12] && t_1<2.0){
	  TF1* parasigT= new TF1("parasigT","expo");
	  parasigT->FixParameter(0, SigPar[1][x_2][y_2][8]); 
	  parasigT->FixParameter(1, SigPar[1][x_2][y_2][9]);
	  sigT14=parasigT->Eval(t_4); 
	  
	  if (SigPar[1][x_2][y_2][8] == 0 && SigPar[1][x_2][y_2][9] == 0){
	    sigT14 =0;
	  }
	}
	 
	else {
	  sigT14=0;
	}

	// Taking the log of claculated sigT values
	lsigT1 = TMath::Log(sigT11); //log value of sigT11.
	lsigT3 = TMath::Log(sigT3); //log value of sigT3.
	lsigT4 = TMath::Log(sigT14); //log value of sigT4.
	// Calculated slopes of different lines
	slb = lsigT4-lsigT3; //->interpolation from the third corner
	sll = -(lsigT1-lsigT3);
      
	// Applying taylor's series formula without averaging the slopes        
	fsigTLa = lsigT3 + (q2-q2_1)*slb + (w-w_2)*sll; //->interpolation from the third corner
  
	// Find the anti-log of the taylor's series formula value 
        fsigTa = exp(fsigTLa);
      }
    //...................................................................................

    else if (sigT1 == 0 && sigT2 == 0){  // if we loose the 1st and the 2nd corner simultaneously
      // In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) and the fourth corner (i.e. top right)will always be there and for other one corner, find the value of the cross-section at the first corner at the minimum value of t. After that we can interpolate them.

      // First try to find t_1

      if(SigPar[1][x_1][y_1][10] != -10001){ // t_1 is t_min for the first corner - this is always the bottom left corner
	t_1 = SigPar[1][x_1][y_1][10];
      }
      else if (SigPar[1][x_1][y_1][11] != -10001){
	t_1 =  SigPar[1][x_1][y_1][11];
      } 
      else if ( SigPar[1][x_1][y_1][12] != -10001){
	t_1 =  SigPar[1][x_1][y_1][12];
      }
      else {
	return -100;
      } 
     
      //Calculating the sigT11 at bottom left corner of the square
      if (t_1>=SigPar[1][x_1][y_1][10]&& t_1<SigPar[1][x_1][y_1][11]){ 
	TF1* parasigT= new TF1("parasigT","pol2");
	parasigT->FixParameter(0, SigPar[1][x_1][y_1][2]); 
	parasigT->FixParameter(1, SigPar[1][x_1][y_1][3]);
	parasigT->FixParameter(2, SigPar[1][x_1][y_1][4]);
	sigT11=parasigT->Eval(t_1);
      }
	 
      else if (t_1>=SigPar[1][x_1][y_1][11] && t_1<SigPar[1][x_1][y_1][12]){
	TF1* parasigT= new TF1("parasigT","pol2");
	parasigT->FixParameter(0,SigPar[1][x_1][y_1][5]); 
	parasigT->FixParameter(1,SigPar[1][x_1][y_1][6]);
	parasigT->FixParameter(2,SigPar[1][x_1][y_1][7]);
	sigT11=parasigT->Eval(t_1); 
      } 
	
      else if (t_1>=SigPar[1][x_1][y_1][12] && t_1<2.0){
	TF1* parasigT= new TF1("parasigT","expo");
	parasigT->FixParameter(0, SigPar[1][x_1][y_1][8]); 
	parasigT->FixParameter(1, SigPar[1][x_1][y_1][9]);
	sigT11=parasigT->Eval(t_1); 	  
	if (SigPar[1][x_1][y_1][8] == 0 && SigPar[1][x_1][y_1][9] == 0){
	  sigT11 =0;
	}
      }
	 
      else {
	sigT11=0;
      }
      // Taking the log of claculated sigT values
       
      lsigT1 = TMath::Log(sigT11); //log value of sigT11.
      lsigT3 = TMath::Log(sigT3); //log value of sigT3.
      lsigT4 = TMath::Log(sigT4); //log value of sigT4.
      // Calculated slopes of different lines
      slb = lsigT3-lsigT4; //->interpolation from the third corner
      sll = -(lsigT1-lsigT3);
      // Applying taylor's series formula without averaging the slopes
      fsigTLa = lsigT3 + (q2-q2_1)*slb + (w-w_2)*sll; //->interpolation from the third corner
        
      // Find the anti-log of the taylor's series formula value         
      fsigTa = exp(fsigTLa);
    }          

    //..............................................................................

    else if (sigT2 == 0) { // if we loose 2nd corner, first we will always loose this corner as this correspond to highest -t value (interpolate from 3rd corner)
      // In this case, we will need atleat three corners to find the cross-section. And even after loosing second corner, we still have three corners to interpolate.

      // Taking the log of claculated sigT values         
      lsigT1 = TMath::Log(sigT1); //log value of sigT1.
      lsigT3 = TMath::Log(sigT3); //log value of sigT3.
      lsigT4 = TMath::Log(sigT4); //log value of sigT4.
    
      // Calculated slopes of different lines
      slb = lsigT3-lsigT4; //->interpolation from the third corner
      sll = -(lsigT1-lsigT3);
        
      // Applying taylor's series formula without averaging the slopes  
      fsigTLa = lsigT3 + (q2-q2_1)*slb + (w-w_2)*sll; //->interpolation from the third corner
      
      // Find the anti-log of the taylor's series formula value       
      fsigTa = exp(fsigTLa);
    }                 
  } // end of if statement over t
  //....................................................................................................................................................................
  else{
    //cerr<<" Invalid t-value "<<endl;
    return -100;
  }      

  //....................................................................................................................................................................
  //.................................................................................................................................................................... 
  // Calcualtion of t-min values at each w and q2 (i.e. at all the four combinations of w's and q2's ) [SigL Calculations started from here]
  //....................................................................................................................................................................
  //.................................................................................................................................................................... 
  // t1 is the t_min value UNLESS the first fit fails
  // If the first fit fails, t2 is the t_min value. If the first fit did NOT fail, t2 is the intersection point between fits 1 and 2
  // If the second fit fails, t3 is the t_min value. If the second fit did NOT fail, t3 is the intersection point between fits 2 and 3
 
  //....................................................................................................................................................................     
  // Variables that need to be defined globally somewhere ( These variables are only for sigL)
  double l_1,l_2,l_3,l_4; // 4 corners (points) in t to interpolate over
  double sigL1,sigL2,sigL3,sigL4; // Value of SigL at each corner of the square to interpolate over
  double lsigL1,lsigL2,lsigL3,lsigL4, stb,stt,stl,str,attb,atlr,fsigLLa,fsigLa; // Logarithm of sigma values and slopes along the bottom, top, left, right and the average of the slope along the top/bottom and left/right. fsigTLa is the result in logarithmic form and fsigTa is the result
  double sigL11,sigL14; // Bottom left and top right sigma points, these are calculated if there are less than 3 (or 2) points to interpolate over
  //...................................................................................................................................................................
       

  /*if(SigPar[0][x_1][y_1][9] != -10001){ // l_1 is t_min for the first corner - this is always the bottom left corner
    l_1 = SigPar[0][x_1][y_1][9];
    }
    else if (SigPar[0][x_1][y_1][10] != -10001){
    l_1 =  SigPar[0][x_1][y_1][10];
    }       
    else if ( SigPar[0][x_1][y_1][11] != -10001){
    l_1 =  SigPar[0][x_1][y_1][11];
    }
    else {
    return -100;
    } 
    //................................................................................... 
  
    if(SigPar[0][x_1][y_2][9] != -10001) {  // l_2 is t_min for the second corner - this is always the bottom right corner
    l_2 = SigPar[0][x_1][y_2][9]; 
    }
    else if (SigPar[0][x_1][y_2][10] != -10001){
    l_2 = SigPar[0][x_1][y_2][10];
    }
    else if (SigPar[0][x_1][y_2][11] != -10001){
    l_2 = SigPar[0][x_1][y_2][11];
    }*/
  //...................................................................................        
 
  if(SigPar[0][x_2][y_1][9] != -10001){ // l_3 is t_min for the third corner - this is always the top left corner
    l_3 = SigPar[0][x_2][y_1][9];
  }
  else if (SigPar[0][x_2][y_1][10] != -10001){
    l_3 = SigPar[0][x_2][y_1][10];
  }
  else if (SigPar[0][x_2][y_1][11] != -10001){
    l_3 = SigPar[0][x_2][y_1][11];
  }
  else {
    return -100;
  } 
  
  //...................................................................................
  
  /*if(SigPar[0][x_2][y_2][9] != -10001){  // l_4 is t_min for the fourth corner - this is always the top right corner
    l_4 = SigPar[0][x_2][y_2][9];
    }
    else if (SigPar[0][x_2][y_2][10] != -10001){
    l_4 =SigPar[0][x_2][y_2][10];
    }
    else if (SigPar[0][x_2][y_2][11] != -10001){
    l_4 = SigPar[0][x_2][y_2][11];
    } */
  //.................................................................................................................................................................... 
  // Calcualtion of sigL's at all the corners of the square 
  //....................................................................................................................................................................  
 
  if (t>= l_3 && t<2.0){ // l_3 corner is the one that we will loose at the end

    //Calculating the sigL1 at bottom left corner of the square..................................................................................
    if (t>=SigPar[0][x_1][y_1][9] && t<SigPar[0][x_1][y_1][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	sigL1=parasigL->Eval(t);
	if (SigPar[0][x_1][y_1][2] == 0 && SigPar[0][x_1][y_1][3] == 0){
	  sigL1=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	parasigL->FixParameter(2, SigPar[0][x_1][y_1][4]);
	sigL1=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_1][10] && t<SigPar[0][x_1][y_1][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_1][5]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_1][6]);
      sigL1=parasigL->Eval(t);
      if (SigPar[0][x_1][y_1][5] == 0 && SigPar[0][x_1][y_1][6] == 0){
	sigL1=0;
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_1][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_1][7]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_1][8]);
      sigL1=parasigL->Eval(t);
      if (SigPar[0][x_1][y_1][7] == 0 && SigPar[0][x_1][y_1][8] == 0){
	sigL1=0;
      }
    }
	 
    else {
      sigL1=0;
    }

    //Calculating the sigL2 at bottom right corner of the square..................................................................................

    if (t>=SigPar[0][x_1][y_2][9] && t<SigPar[0][x_1][y_2][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_2][3]);
	sigL2=parasigL->Eval(t);
	if (SigPar[0][x_1][y_2][2] == 0 && SigPar[0][x_1][y_2][3] == 0){
	  sigL2=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_1][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_2][3]);
	parasigL->FixParameter(2, SigPar[0][x_1][y_2][4]);
	sigL2=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_2][10] && t<SigPar[0][x_1][y_2][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_2][5]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_2][6]);
      sigL2=parasigL->Eval(t); 
      if (SigPar[0][x_1][y_2][5] == 0 && SigPar[0][x_1][y_2][6] == 0){
	sigL2=0;
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_2][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_2][7]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_2][8]);
      sigL2=parasigL->Eval(t);
      if (SigPar[0][x_1][y_2][7] == 0 && SigPar[0][x_1][y_2][8] == 0){
	sigL2=0;
      }
    }
	 
    else {
      sigL2=0;
    }

    //Calculating the sigL3 at top left corner of the square..................................................................................
    if (t>=SigPar[0][x_2][y_1][9] && t<SigPar[0][x_2][y_1][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_2][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_1][3]);
	sigL3=parasigL->Eval(t); 
	 
	if (SigPar[0][x_2][y_1][2] == 0 && SigPar[0][x_2][y_1][3] == 0){
	  sigL3=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_2][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_1][3]);
	parasigL->FixParameter(2, SigPar[0][x_2][y_1][4]);
	sigL3=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_1][10] && t<SigPar[0][x_2][y_1][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_1][5]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_1][6]);
      sigL3=parasigL->Eval(t); 
      
      if (SigPar[0][x_2][y_1][5] == 0 && SigPar[0][x_2][y_1][6] == 0){
	sigL3=0;
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_1][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_1][7]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_1][8]);
      sigL3=parasigL->Eval(t); 
	 
      if (SigPar[0][x_2][y_1][7] == 0 && SigPar[0][x_2][y_1][8] == 0){
	sigL3=0;
      }
    }
	 
    else {
      sigL3=0;
    }

    //Calculating the sigL4 at top right corner of the square..................................................................................
    if (t>=SigPar[0][x_2][y_2][9] && t<SigPar[0][x_2][y_2][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	sigL4=parasigL->Eval(t); 
	 
	if (SigPar[0][x_2][y_2][2] == 0 && SigPar[0][x_2][y_2][3] == 0){
	  sigL4=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	parasigL->FixParameter(2, SigPar[0][x_2][y_2][4]);
	sigL4=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_2][10] && t<SigPar[0][x_2][y_2][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_2][5]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_2][6]);
      sigL4=parasigL->Eval(t); 
      
      if (SigPar[0][x_2][y_2][5] == 0 && SigPar[0][x_2][y_2][6] == 0){
	sigL4=0;
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_2][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_2][7]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_2][8]);
      sigL4=parasigL->Eval(t); 
	 
      if (SigPar[0][x_2][y_2][7] == 0 && SigPar[0][x_2][y_2][8] == 0){
	sigL4=0;
      }
    }
	 
    else {
      sigL4=0;
    }

    //....................................................................................................................................................................
    // Difterent if and else conditions to find the crossection values at the given points.
    //....................................................................................................................................................................

    if (sigL1 == sigL2 && sigL2 == sigL3 && sigL3 == sigL4 && sigL4 == sigL1){ // if the w and q2 will have whole number values
      cerr<<"fsigLa = "<<sigL1 <<endl;
    }   

    //...................................................................................

    else if (sigL1 != 0 && sigL2 != 0 && sigL3 != 0 && sigL4 != 0){ // if all the four corners are  present
    
      // Taking the log of claculated sigT values
      lsigL1 = TMath::Log(sigL1); //log value of sigT1.
      lsigL2 = TMath::Log(sigL2); //log value of sigT2.
      lsigL3 = TMath::Log(sigL3); //log value of sigT3.
      lsigL4 = TMath::Log(sigL4); //log value of sigT4.
     
      // Calculated slopes of different lines
       
      stb = lsigL4-lsigL3; //interpolation from the third corner (i.e. top left corner)
      stt = lsigL2-lsigL1;
      stl = -(lsigL1-lsigL3);
      str = -(lsigL2-lsigL4);
        
      // Taking averages of the slopes
        
      attb = ((stb +stt)/2 );
      atlr = ((stl +str)/2 );
        
      //  Applying taylor's series formula for the average slopes
        
      fsigLLa = lsigL3 + (q2-q2_1)*attb + (w-w_2)*atlr;
        
      // Find the anti-log of the taylor's series formula value
      fsigLa = exp(fsigLLa);
    }
    //...................................................................................

    else if (sigL1 == 0 && sigL2 == 0 && sigL4 == 0)  // if we loose the 1st, 2nd and the 4th corner simultaneously 
      {
	// In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) will always be there and for other two corners, find the value of the cross-section at the first and the fourth corner at the minimum value of t. After that we can interpolate them.
 
	// First try to find t_1 and t_4

	if(SigPar[0][x_1][y_1][9] != -10001){ // l_1 is t_min for the first corner - this is always the bottom left corner
	  l_1 = SigPar[0][x_1][y_1][9];
	}
	else if (SigPar[0][x_1][y_1][10] != -10001){
	  l_1 =  SigPar[0][x_1][y_1][10];
	}       
	else if ( SigPar[0][x_1][y_1][11] != -10001){
	  l_1 =  SigPar[0][x_1][y_1][11];
	}
	else {
	  return -100;
	} 
	//...................................................................................
	if(SigPar[0][x_2][y_2][9] != -10001){  // l_4 is t_min for the fourth corner - this is always the top right corner
	  l_4 = SigPar[0][x_2][y_2][9];
	}
	else if (SigPar[0][x_2][y_2][10] != -10001){
	  l_4 =SigPar[0][x_2][y_2][10];
	}
	else if (SigPar[0][x_2][y_2][11] != -10001){
	  l_4 = SigPar[0][x_2][y_2][11];
	}
	else {
	  return -100;
	} 
	//Calculating the sigL11 at bottom left corner of the square
	if (l_1>=SigPar[0][x_1][y_1][9] && l_1<SigPar[0][x_1][y_1][10]){
	  if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	    TF1* parasigL= new TF1("parasigL","expo");
	    parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	    sigL11=parasigL->Eval(l_1); 
	 
	    if (SigPar[0][x_1][y_1][2] == 0 && SigPar[0][x_1][y_1][3] == 0){
	      sigL11=0;
	    }
	  }
      
	  else{
	    TF1* parasigL= new TF1("parasigL","pol2");
	    parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	    parasigL->FixParameter(2, SigPar[0][x_1][y_1][4]);
	    sigL11=parasigL->Eval(l_1);
	  }
	}
	 
	else if (l_1>=SigPar[0][x_1][y_1][10] && l_1<SigPar[0][x_1][y_1][11]){
	  TF1* parasigL= new TF1("parasigL","expo");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][5]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][6]);
	  sigL11=parasigL->Eval(l_1); 
      
	  if (SigPar[0][x_1][y_1][5] == 0 && SigPar[0][x_1][y_1][6] == 0){
	    sigL11=0;
	  }
	}
	 
	else if (l_1>=SigPar[0][x_1][y_1][11] && t<2.0){
	  TF1* parasigL= new TF1("parasigT","expo");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][7]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][8]);
	  sigL11=parasigL->Eval(l_1); 
	 
	  if (SigPar[0][x_1][y_1][7] == 0 && SigPar[0][x_1][y_1][8] == 0){
	    sigL11=0;
	  }
	}
	 
	else {
	  sigL11=0;
	}   
    
	//Calculating the sigL14 at bottom left corner of the square 
	if (l_4>=SigPar[0][x_2][y_2][9] && l_4<SigPar[0][x_2][y_2][10]){
	  if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	    TF1* parasigL= new TF1("parasigL","expo");
	    parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	    sigL14=parasigL->Eval(l_4); 
	 
	    if (SigPar[0][x_2][y_2][2] == 0 && SigPar[0][x_2][y_2][3] == 0){
	      sigL14=0;
	    }
	  }
      
	  else{
	    TF1* parasigL= new TF1("parasigL","pol2");
	    parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	    parasigL->FixParameter(2, SigPar[0][x_2][y_2][4]);
	    sigL14=parasigL->Eval(l_4);
	  }
	}
	 
	else if (l_4>=SigPar[0][x_2][y_2][10] && l_4<SigPar[0][x_2][y_2][11]){
	  TF1* parasigL= new TF1("parasigL","expo");
	  parasigL->FixParameter(0, SigPar[0][x_2][y_2][5]); 
	  parasigL->FixParameter(1, SigPar[0][x_2][y_2][6]);
	  sigL14=parasigL->Eval(l_4); 
      
	  if (SigPar[0][x_2][y_2][5] == 0 && SigPar[0][x_2][y_2][6] == 0){
	    sigL14=0;
	  }
	}
	 
	else if (l_4>=SigPar[0][x_2][y_2][11] && t<2.0){
	  TF1* parasigL= new TF1("parasigT","expo");
	  parasigL->FixParameter(0, SigPar[0][x_2][y_2][7]); 
	  parasigL->FixParameter(1, SigPar[0][x_2][y_2][8]);
	  sigL14=parasigL->Eval(l_4); 
	 
	  if (SigPar[0][x_2][y_2][7] == 0 && SigPar[0][x_2][y_2][8] == 0){
	    sigL14=0;
	  }
	}
	 
	else {
	  sigL14=0;
	}

	// Taking the log of claculated sigL values
	lsigL1 = TMath::Log(sigL11); //log value of sigL11.
	lsigL3 = TMath::Log(sigL3); //log value of sigL3.
	lsigL4 = TMath::Log(sigL14); //log value of sigL14.

	// Calculated slopes of different lines  
	stb = lsigL4-lsigL3; //->interpolation from the third corner
	stl = -(lsigL1-lsigL3);
	// Applying taylor's series formula without averaging the slopes
     	fsigLLa = lsigL3 + (q2-q2_1)*stb + (w-w_2)*stl; //->interpolation from the third corner
     
	// Find the anti-log of the taylor's series formula value 
        fsigLa = exp(fsigLLa);
      }
    //...................................................................................

    else if (sigL1 == 0 && sigL2 == 0){  // if we loose the 1st and the 2nd corner simultaneously
      // In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) and the fourth corner (i.e. top right)will always be there and for other one corner, find the value of the cross-section at the first corner at the minimum value of t. After that we can interpolate them.

      // First try to find t_1

      if(SigPar[0][x_1][y_1][9] != -10001){ // l_1 is t_min for the first corner - this is always the bottom left corner
	l_1 = SigPar[0][x_1][y_1][9];
      }
      else if (SigPar[0][x_1][y_1][10] != -10001){
	l_1 =  SigPar[0][x_1][y_1][10];
      }       
      else if ( SigPar[0][x_1][y_1][11] != -10001){
	l_1 =  SigPar[0][x_1][y_1][11];
      }
      else {
	return -100;
      } 
      //...................................................................................
      //Calculating the sigL11 at bottom left corner of the square
      if (l_1>=SigPar[0][x_1][y_1][9] && l_1<SigPar[0][x_1][y_1][10]){
	if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	  TF1* parasigL= new TF1("parasigL","expo");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	  sigL11=parasigL->Eval(l_1); 
	 
	  if (SigPar[0][x_1][y_1][2] == 0 && SigPar[0][x_1][y_1][3] == 0){
	    sigL11=0;
	  }
	}
      
	else{
	  TF1* parasigL= new TF1("parasigL","pol2");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	  parasigL->FixParameter(2, SigPar[0][x_1][y_1][4]);
	  sigL11=parasigL->Eval(l_1);
	}
      }
	 
      else if (l_1>=SigPar[0][x_1][y_1][10] && l_1<SigPar[0][x_1][y_1][11]){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][5]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][6]);
	sigL11=parasigL->Eval(l_1);
	if (SigPar[0][x_1][y_1][5] == 0 && SigPar[0][x_1][y_1][6] == 0){
	  sigL11=0;
	}
      }
	 
      else if (l_1>=SigPar[0][x_1][y_1][11] && t<2.0){
	TF1* parasigL= new TF1("parasigT","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][7]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][8]);
	sigL11=parasigL->Eval(l_1);	 
	if (SigPar[0][x_1][y_1][7] == 0 && SigPar[0][x_1][y_1][8] == 0){
	  sigL11=0;
	}
      }
	 
      else {
	sigL11=0;
      }
     
      // Taking the log of claculated sigL values
      lsigL1 = TMath::Log(sigL11); //log value of sigL11.
      lsigL3 = TMath::Log(sigL3); //log value of sigL3.
      lsigL4 = TMath::Log(sigL4); //log value of sigL4.
      // Calculated slopes of different lines 
      stb = lsigL3-lsigL4; //->interpolation from the third corner
      stl = -(lsigL1-lsigL3);
        
      // Applying taylor's series formula without averaging the slopes
      fsigLLa = lsigL3 + (q2-q2_1)*stb + (w-w_2)*stl; //->interpolation from the third corner
      
      // Find the anti-log of the taylor's series formula value        
      fsigLa = exp(fsigLLa);
    }          

    //..............................................................................

    else if (sigL2 == 0) { // if we loose 2nd corner, first we will always loose this corner as this correspond to highest -t value (interpolate from 3rd corner)
      // In this case, we will need atleat three corners to find the cross-section. And even after loosing second corner, we still have three corners to interpolate.      
      // Taking the log of claculated sigT values         
      lsigL1 = TMath::Log(sigL1); //log value of sigL1.
      lsigL3 = TMath::Log(sigL3); //log value of sigL3.
      lsigL4 = TMath::Log(sigL4); //log value of sigL4.

      // Calculated slopes of different lines
      stb = lsigL3-lsigL4; //->interpolation from the third corner
      stl = -(lsigL1-lsigL3);
        
      // Applying taylor's series formula without averaging the slopes
      fsigLLa = lsigL3 + (q2-q2_1)*stb + (w-w_2)*stl; //->interpolation from the third corner
   
      // Find the anti-log of the taylor's series formula value         
      fsigLa = exp(fsigLLa);
    }                 
  } // end of if statement over t
  //....................................................................................................................................................................
  else{
    //cerr<<" Invalid t-value "<<endl;
    return -100;
  } 
  //....................................................................................................................................................................
  //.................................................................................................................................................................... 

  //sig_total = 5222;
  sig_total = fsigTa +(feps*fsigLa); 

  return sig_total;

}

// SJDK 21/12/22 - This function needs updating!
Double_t DEMP_Reaction::GetPi0_CrossSection() {

  double_t sig_total;
  return sig_total;

}

/*--------------------------------------------------*/
/// Output generator detail

void DEMP_Reaction::Detail_Output() {

  DEMPDetails << "Total events tried                                           " << setw(20) << fNGenerated   << endl;
  DEMPDetails << "Total events recorded                                        " << setw(20) << fNRecorded    << endl;
  DEMPDetails << "Number of events with wsq negative                           " << setw(20) << w_neg_ev      << endl;
  DEMPDetails << "Number of events with 2 < w < 10                             " << setw(20) << w_ev          << endl;
  DEMPDetails << "Number of events with qsq < 1 or > 35                        " << setw(20) << qsq_ev        << endl;
  DEMPDetails << "Number of events with Meson (X) energy NaN                   " << setw(20) << fNaN          << endl;
  DEMPDetails << "Number of events failing conservation law check              " << setw(20) << fConserve     << endl;
  DEMPDetails << "Total events passing conservation laws                       " << setw(20) << conserve   << endl;
  DEMPDetails << "Total events failed energy conservation                      " << setw(20) << ene   << endl; 
  DEMPDetails << "Total events failed momentum conserveation                   " << setw(20) << mom   << endl;
  DEMPDetails << "Number of events with -t > 2 (K+) or -t > 1.3 (Pi+) GeV      " << setw(20) << t_ev          << endl;
  DEMPDetails << "Number of events with w less than threshold                  " << setw(20) << fWSqNeg       << endl;
  DEMPDetails << "Number of events with mom not conserve                       " << setw(20) << fNMomConserve << endl;
  DEMPDetails << "Number of events with Sigma negative                         " << setw(20) << fNSigmaNeg    << endl;
  DEMPDetails << "Number of lund events                                        " << setw(20) << fLundRecorded << endl;

  DEMPDetails << "Seed used for the Random Number Generator                    " << setw(20) << fSeed         << endl;

}

////*--------------------------------------------------
/// Functions for different output formats follow

void DEMP_Reaction::Lund_Output() {

  DEMPOut << "3"
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
  DEMPOut << setw(10) << "1" 
	  << setw(10) << "1" 
	  << setw(10) << "1" 
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
  DEMPOut << setw(10) << "2" 
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
 	  
  // Recoiled neutron
  DEMPOut << setw(10) << "3" 
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

void DEMP_Reaction::DEMPReact_Pythia6_Out_Init() {

	print_itt = 0;

	DEMPOut << "DEMP Event FILE" << endl;
	DEMPOut << "============================================" << endl;
	DEMPOut << "I, ievent, nParticles, Weight" << endl;
	DEMPOut << "============================================" << endl;
	DEMPOut << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)" << endl;
	DEMPOut << "============================================" << endl;

}

void DEMP_Reaction::DEMPReact_Pythia6_Output() {

  DEMPOut << "0" << " \t\t\t "  << print_itt << " \t\t\t " << "1" << " \t\t\t " << fEventWeight << endl; // var 1

  print_itt++;

  DEMPOut << "============================================" << endl;

  ///*--------------------------------------------------*/
  // Initial State
 
  DEMPOut  << "1" 
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

  DEMPOut << "2" 
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

  DEMPOut << "3" 
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
  DEMPOut << "4" 
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
  DEMPOut << "5" 
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
  DEMPOut << "6" 
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

  DEMPOut << "=============== Event finished ===============" << endl;

}

/*--------------------------------------------------*/

void DEMP_Reaction::DEMPReact_HEPMC3_Out_Init() {
 
  print_itt = 0;
  DEMPOut << "HepMC::Version 3.02.02" << endl;
  DEMPOut << "HepMC::Asciiv3-START_EVENT_LISTING" << endl;

}

/*--------------------------------------------------*/

void DEMP_Reaction::DEMPReact_HEPMC3_Output() {
  
  // HEPMC3 output for Athena/EPIC simulations

  // First line - E - Event# - #Vertices - #Particles
  DEMPOut << "E" << " "  << print_itt <<  " " << "1" << " " << 5 << endl;
  print_itt++;
  // Second line, Units - U - ENERGY UNIT - DISTANCE UNIT
  DEMPOut << "U" << " " << "GEV" << " " << "MM" << endl;
  // Third line, optional attributes, the weight
  DEMPOut << "A" << " " << "0" << " " << "weight" << " " <<  fEventWeight << endl;
  // Beam particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (4, incoming beam particle)
  DEMPOut << "P" << " " << "1" << " " << "0" << " " << "11" << " " << r_lelectrong.X() << " " << r_lelectrong.Y() << " " << r_lelectrong.Z() << " " << r_lelectrong.E() << " " << fElectron_Mass_GeV << " " << "4" << endl;
  DEMPOut << "P" << " " << "2" << " " << "0" << " " << "2212" << " " << r_lprotong.X() << " " << r_lprotong.Y() << " " << r_lprotong.Z() << " " << r_lprotong.E() << " " << fProton_Mass_GeV << " " << "4" << endl;
  // Vertex line - V - 1 - 0 - [1,2]
  DEMPOut << "V" << " " << "-1" << " " << "0" << " " << "[1,2]" << endl;
  // Output particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (1, undecayed physical particle)
  // Scattered electron
  DEMPOut << "P" << " " << "3" << " " << "-1" << " " << "11" << " " << r_lscatelecg.X() << " "  << r_lscatelecg.Y() << " "  << r_lscatelecg.Z() << " " << r_lscatelecg.E() << " " << fElectron_Mass_GeV << " " << "1" << endl;
  // Produced meson
  DEMPOut << "P" << " " << "4" << " " << "-1" << " " << PDGtype(produced_X) << " " << r_lX_g.X() << " "  << r_lX_g.Y() << " "  << r_lX_g.Z() << " " << r_lX_g.E() << " " << fX_Mass_GeV << " " << "1" << endl;
  // Recoil hadron
  DEMPOut << "P" << " " << "5" << " " << "-1" << " " << PDGtype(recoil_hadron) << " " << r_l_scat_hadron_g.X() << " "  << r_l_scat_hadron_g.Y() << " "  << r_l_scat_hadron_g.Z() << " " << r_l_scat_hadron_g.E() << " " << f_Scat_hadron_Mass_GeV << " " << "1" << endl;
  
}
