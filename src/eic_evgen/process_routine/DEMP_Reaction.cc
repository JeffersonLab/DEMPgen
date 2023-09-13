#include "reaction_routine.h"
#include "eic.h"

#include <sys/stat.h>

using namespace std;

DEMP_Reaction::DEMP_Reaction() { 

  cout << "Program Start" << endl;

}

/*--------------------------------------------------*/
/// DEMP_Reaction 

DEMP_Reaction::DEMP_Reaction(TString particle_str, TString hadron_str) { 

  rEjectile = particle_str;
  rRecoil = hadron_str;

}

DEMP_Reaction::~DEMP_Reaction() {

  // Text output data file and generation information 

  DEMPOut.close();
  DEMPDetails.close();

  // Diagnostic root file with plots

  dRootTree->Write();
  dRootFile->Close();

  delete dRootFile;
  delete dRootTree;

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
 
    Progress_Report();  // This happens at each 10% of the total event is processed
    Processing_Event();
  }
 
  Detail_Output();
 
}

void DEMP_Reaction::Init() {

  pim* myPim;

  pd = dynamic_cast<pim*>(myPim);
	
  rEjectile_charge = ExtractCharge(rEjectile);

  const char* dir_name = "OutputFiles";
  struct stat sb;

  if (stat(dir_name, &sb) == 0) {
     cout << "The path is valid!";
  }
  else {
     cout << "The Path is invalid!";
     mkdir(dir_name,0777);
  } 

  sTFile = Form("./%s/eic_%s.txt", dir_name, gfile_name.Data());
  sLFile = Form("./%s/eic_input_%s.dat", dir_name, gfile_name.Data());
  sDFile = Form("./%s/eic_input_%s.root", dir_name, gfile_name.Data());

  DEMPOut.open( sLFile.c_str() );
  DEMPDetails.open( sTFile.c_str() );
 
  dRootFile = new TFile(sDFile.c_str(),"RECREATE"); 
  dRootTree = new TTree();
 
  dRootTree->Branch("EventWeight", &fEventWeight, "fEventWeight/D");

  /*--------------------------------------------------*/ 
	
  qsq_ev = 0, t_ev = 0, w_neg_ev = 0, w_ev = 0;
  rNEvents = fNEvents;
  rNEvent_itt = 0;
  
  // 02/06/21 SJDK 
  // Set these values once the beam energies are read in
  fPSF = ( fEBeam * ( fScatElec_E_Hi - fScatElec_E_Lo ) *( sin( fScatElec_Theta_F ) - sin( fScatElec_Theta_I ) ) * 2 * fPI *( sin( fEjectileX_Theta_F ) - sin( fEjectileX_Theta_I ) ) * 2 * fPI );
  fElectron_Kin_Col_GeV = fEBeam;
  fElectron_Kin_Col = fElectron_Kin_Col_GeV * 1000.0;

  // cout << rNEvents << "    " << fNEvents << endl;

  // 08/09/23 - SJDK - Fermi momentum commented out for now, this is not fully implemented yet
  // In future, this will be enabled/disabled automatically depending upon the specified hadron beam
  //rFermiMomentum = pd->fermiMomentum();

  // ----------------------------------------------------
  // Proton in collider (lab) frame
  // 08/09/23 - SJDK - The naming needs to be adjusted to be the incoming hadron beam, not the proton. Make this generic.
  r_lproton = GetProtonVector_lab();
  r_lprotong = GetProtonVector_lab() * fm;

  // Getting the mass of the hadron beam
  r_lhadron_beam_mass = ParticleMass(ParticleEnum(gBeamPart.c_str()))*1000; // in MeV

//  cout << gBeamPart << endl;
//  cout << ParticleEnum(gBeamPart.c_str()) << endl; 
//  cout << ParticleMass(ParticleEnum(gBeamPart.c_str())) << endl; 
//  cout << r_lhadron_beam_mass << endl;
//  exit(0);

  // ----------------------------------------------------
  // Electron in collider (lab) frame

  // 06/09/23 - SJDK - Commenting out for now, should be disabled for e/p collisions
  //cout << "Fermi momentum: " << rFermiMomentum << endl;

  r_lelectron	 = GetElectronVector_lab();
  r_lelectrong = r_lelectron * fm;

  ///*--------------------------------------------------*/
  /// Getting the ejectile (produced meson) particle mass from the data base
 
  produced_X = ParticleEnum(rEjectile);
  f_Ejectile_Mass = ParticleMass(produced_X)*1000; //MeV
  f_Ejectile_Mass_GeV = f_Ejectile_Mass/1000; //GeV

  cout << rEjectile << "  " << produced_X << "  " << f_Ejectile_Mass_GeV <<  endl;
  cout << rEjectile_charge << endl;

  if (rRecoil == "Neutron" ) {
    rEjectile_scat_hadron  = "Neutron"; 
    recoil_hadron  = Neutron; 
    f_Recoil_Mass     = fNeutron_Mass;
    f_Recoil_Mass_GeV = f_Recoil_Mass/1000;
  }
  else if (rRecoil == "Proton" ) {	
    rEjectile_scat_hadron  = "Proton"; 
    recoil_hadron  = Proton;
    f_Recoil_Mass = fProton_Mass;
    f_Recoil_Mass_GeV = f_Recoil_Mass/1000;
  } 
  else if(rRecoil == "Lambda"){
    rEjectile_scat_hadron  = "Lambda"; 
    recoil_hadron  = Lambda;
    f_Recoil_Mass = fLambda_Mass;
    f_Recoil_Mass_GeV = f_Recoil_Mass/1000;
    cout<<"Particle = "<<rEjectile_scat_hadron<<" with mass = "<< f_Recoil_Mass << endl;
  }
  else if (rRecoil == "Sigma0"){
    rEjectile_scat_hadron  = "Sigma0";
    recoil_hadron  = Sigma0;
    f_Recoil_Mass  = fSigma_Mass;
    f_Recoil_Mass_GeV = f_Recoil_Mass/1000;
    cout<<"Particle = "<<rEjectile_scat_hadron<<" with mass = "<< f_Recoil_Mass << endl;
  }

  rDEG2RAD   = fPI/180.0;

  f_Ejectile_Theta_I = fEjectileX_Theta_I ;
  f_Ejectile_Theta_F = fEjectileX_Theta_F;

  cout << "Produced particle in exclusive production: " << rEjectile << ";  with mass: " << f_Ejectile_Mass << " MeV "<< endl;
  cout << fEBeam << " GeV electrons on " << fHBeam << " GeV ions" << endl;
  if(UseSolve == true){
    cout << rEjectile << " and " << rEjectile_scat_hadron << " 4-vectors calculated using Solve function" << endl;
  }
  else if(UseSolve == false){
    cout << rEjectile << " and " << rEjectile_scat_hadron << " 4-vectors calculated using analytical solution" << endl;
  }
  // Set luminosity value based upon beam energy combination, note that if no case matches, a default of 1e33 is assumed. Cases are a set of nominal planned beam energy combinations for the EIC (and EICC)
  // See slide 11 in https://indico.cern.ch/event/1072579/contributions/4796856/attachments/2456676/4210776/CAP-EIC-June-7-2022-Seryi-r2.pdf
  // If available in the future, this could be replaced by some fixed function
  if ((fEBeam == 5.0 ) && (fHBeam == 41.0) ){
    fLumi = 0.44e33;
  }
  else if ((fEBeam == 5.0 ) && (fHBeam == 100.0) ){
    fLumi = 3.68e33;
  }
  else if ((fEBeam == 10.0 ) && (fHBeam == 100.0) ){
    fLumi = 4.48e33;
  }
  else if ((fEBeam == 18.0 ) && (fHBeam == 275.0) ){
    fLumi = 1.54e33;
  }
  else if ((fEBeam == 3.5 ) && (fHBeam == 20) ){ // EICC optimal beam energy combination
    fLumi = 2e33;
  }
  else if ((fEBeam == 2.8 ) && (fHBeam == 13) ){ // EICC lowest beam energy combination
    fLumi = 0.7e33;
  }
  else{
    cout << "!!! Notice !!! The beam energy combination simulated does not match an expected case, a default luminosity value of - " << fLumi << " cm^2s^-1 has been assumed. !!! Notice !!!" << endl;
  }

  if(UseSolve == true){
    /*--------------------------------------------------*/ 
    // SJDK 03/04/22 -  New set of initialisation stuff for the solve function from Ishan and Bill
    
    CoinToss = new TRandom3();

    F = new TF1("F",
		"[6]-sqrt([7]**2+x**2)-sqrt([8]**2+([3]-[0]*x)**2+([4]-[1]*x)**2+([5]-[2]*x)**2)",
		0, r_lproton.E());

    char AngleGenName[100] = "AngleGen";
    double dummy[2] = {0,1};

    f_Ejectile_Theta_I = fEjectileX_Theta_I;
    f_Ejectile_Theta_F = fEjectileX_Theta_F;

    double ThetaRange[2] = {f_Ejectile_Theta_I, f_Ejectile_Theta_F};
    double PhiRange[2] = {0, 360*TMath::DegToRad()};

    AngleGen = new CustomRand(AngleGenName, dummy,
			      ThetaRange, PhiRange);

    UnitVect = new TVector3(0,0,1);

    //  ///*--------------------------------------------------*/ 
    //  // Produced hadron and recoilded hadron from the solve function 

    VertBeamElec = new TLorentzVector();
    VertScatElec = new TLorentzVector();

    Initial      = new TLorentzVector();
    Target       = new TLorentzVector();
    Photon       = new TLorentzVector();
    Interaction  = new TLorentzVector();
    Final        = new TLorentzVector();
  }
  
}

void DEMP_Reaction::Processing_Event() {

  // ----------------------------------------------------
  // Considering Fermi momentum for the proton
  // ----------------------------------------------------
  // SJDK - 06/06/23 - Commenting out for now, this increases the number of recorded events - Should have no fermi momentum for e/p collisions
  // if( kCalcFermi ) {
  //   Consider_Proton_Fermi_Momentum(); 
  // }

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
  // Produced ejectile in Collider frame
  // ----------------------------------------------------  

  f_Ejectile_Theta_Col      = acos( fRandom->Uniform( cos(f_Ejectile_Theta_I), cos(f_Ejectile_Theta_F ) ) ); 
  f_Ejectile_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi );
    	
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
  // SJDK 03/04/23 - Qsq an W ranges now variables set by particle type in .json read in. See eic.cc
  if ( fQsq_GeV < fQsq_Min || fQsq_GeV > fQsq_Max ) {
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

  if ( fW_GeV < fW_Min || fW_GeV > fW_Max ) { // SJDK 03/04/23 - Switched to the new variable, set by particle type
    w_ev++;
    return;
  }

  // SJDK - 06/09/23 - Check UseSolved boolean, process through relevant loop
  if (UseSolve == false){
    // ---------------------------------------------------------
    // Pion momentum in collider frame, analytic solution starts
    // ---------------------------------------------------------
 
    double fupx = sin( f_Ejectile_Theta_Col ) * cos( f_Ejectile_Phi_Col );
    double fupy = sin( f_Ejectile_Theta_Col ) * sin( f_Ejectile_Phi_Col );
    double fupz = cos( f_Ejectile_Theta_Col );
 
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
     
    double ft = fc * fc - fb + f_Ejectile_Mass * f_Ejectile_Mass - f_Recoil_Mass * f_Recoil_Mass;
     
    double fQA = 4.0 * ( fa * fa - fc * fc );
    double fQB = 4.0 * fc * ft;

    double fQC = -4.0 * fa * fa * f_Ejectile_Mass * f_Ejectile_Mass - ft * ft;    
 
    fradical = fQB * fQB - 4.0 * fQA * fQC;
 
    fepi1 = ( -fQB - sqrt( fradical ) ) / ( 2.0 * fQA );
    fepi2 = ( -fQB + sqrt( fradical ) ) / ( 2.0 * fQA );

    ///---------------------------------------------------------
    /// Particle X momentum in collider frame, analytic solution
    /// And obtain recoiled proton in collider (lab) frame
    ///---------------------------------------------------------

    r_l_Ejectile.SetPxPyPzE( (sqrt( pow( fepi1 , 2) - pow(f_Ejectile_Mass , 2) ) ) * sin(f_Ejectile_Theta_Col) * cos(f_Ejectile_Phi_Col),
			     ( sqrt( pow( fepi1 , 2) - pow(f_Ejectile_Mass , 2) ) ) * sin(f_Ejectile_Theta_Col) * sin(f_Ejectile_Phi_Col),
			     ( sqrt( pow( fepi1 , 2) - pow(f_Ejectile_Mass , 2) ) ) * cos(f_Ejectile_Theta_Col),
			     fepi1 );
  
    l_Recoil.SetPxPyPzE( ( r_lproton + r_lelectron - r_lscatelec - r_l_Ejectile).X(),
			 ( r_lproton + r_lelectron - r_lscatelec - r_l_Ejectile ).Y(),
			 ( r_lproton + r_lelectron - r_lscatelec - r_l_Ejectile ).Z(),
			 sqrt( pow( ( ( ( r_lproton + r_lelectron - r_lscatelec - r_l_Ejectile ).Vect() ).Mag()),2) +
			       pow( f_Recoil_Mass , 2) ) );
  }
   else if (UseSolve == true){
     if(!Solve()){
       return;
     }  
     r_l_Ejectile = r_l_Ejectile_solved;
     l_Recoil = r_l_Recoil_solved;
   }
  
  ///--------------------------------------------------
  
  r_l_Ejectile_g = r_l_Ejectile * fm;
  l_Recoil_g = l_Recoil * fm;

  // ----------------------------------------------------------------------------------------------
  // Calculate w = (proton + photon)^2
  // ----------------------------------------------------------------------------------------------
  
  r_lw = r_lproton + r_lphoton;
  fW = r_lw.Mag();

  // SJDK 15/06/21 - Added integer counters for conservation law check and for NaN check
  if (r_l_Ejectile.E() != r_l_Ejectile.E()){ // SJDK 15/06/21 - If the energy of the produced meson is not a number, return and add to counter
    fNaN++;
    return;
  }

  ///*--------------------------------------------------*/ 
  //-> 10/05/23 - Love added a slimmed down, simpler to read version of the CheckLaws fn
  // 
  // To check the conservation of the energy and momentum, there two methods avalaible:
  // Method 1: Give the four-vectors of the initial and final states partciles, 
  //           tolerance factor will be defaulted 1e-6 MeV
  //           CheckLaws(e_beam, h_beam, scatt_e, ejectile, recoil) <- input 4 vectors
  // Method 2: Give the four-vectors of the initial and final states partciles, 
  //           and the prefered tolerance factor.
  //           CheckLaws(e_beam, h_beam, scatt_e, ejectile, recoil, tolerance) <- input 4 vectors and tolerance value in GeV
  // Both functions return 1 if conservations laws are satisified
  
   if( pd->CheckLaws(r_lelectron, r_lproton, r_lscatelec, r_l_Ejectile, l_Recoil) !=1 ){
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
     
  l_Ejectile_rf = r_l_Ejectile;
  l_Ejectile_rf.Boost(-beta_col_rf);
  l_Ejectile_rfg = l_Ejectile_rf * fm;
         
  l_scat_hadron_rf = l_Recoil;
  l_scat_hadron_rf.Boost(-beta_col_rf);
  l_scat_hadron_rf_g = l_scat_hadron_rf * fm;

  ////////////////////////////////////////////////////////////////////////////////////////////
  //                                          End                                           //
  // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
  ////////////////////////////////////////////////////////////////////////////////////////////

  // -----------------------------------------------------------------------------------------
  // Calculate -t
  // -----------------------------------------------------------------------------------------
  // 18/05/23 - SJDK - Should these be proton mass still or not?

  fBeta_CM_RF        = (lphoton_rf.Vect()).Mag() / (lphoton_rf.E() + r_lhadron_beam_mass);

  fGamma_CM_RF       = (lphoton_rf.E() + r_lhadron_beam_mass) / fW;
  f_Ejectile_Energy_CM       = (pow(fW, 2) + pow(f_Ejectile_Mass,2) - pow(f_Recoil_Mass,2) ) / (2.0* fW);    
  f_Ejectile_Mom_CM          = sqrt(pow(f_Ejectile_Energy_CM,2) - pow(f_Ejectile_Mass,2));    
  f_Ejectile_Energy_CM_GeV   = f_Ejectile_Energy_CM / 1000.0;
  f_Ejectile_Mom_CM_GeV      = f_Ejectile_Mom_CM / 1000.0;

  // this equation is valid for parallel kinematics only!
  fT_Para = ( pow(((r_lphoton.Vect()).Mag() - (r_l_Ejectile.Vect()).Mag()),2) - pow((r_lphoton.E() - r_l_Ejectile.E()),2));
  fT_Para_GeV = fT_Para/1000000.0;

  lt = r_lphoton - r_l_Ejectile;
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

  // 31/01/23 SJDK - New limit on t, remove only events outside the parameterisation range
  // 06/09/23 SJDK - fT_Max set in eic.cc depending upon ejectile type
  if (fT_GeV > fT_Max ) {
    t_ev++;
    return;
  }
  
  fx = fQsq_GeV / ( 2.0 * r_lprotong.Dot( r_lphotong ) );
  fy = r_lprotong.Dot( r_lphotong ) / r_lprotong.Dot( r_lelectrong );
  fz = r_l_Ejectile.E()/r_lphoton.E();    

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

  v3X.SetX( l_Ejectile_rfg.X() ) ;        
  v3X.SetY( l_Ejectile_rfg.Y() ) ;        
  v3X.SetZ( l_Ejectile_rfg.Z() );

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
 
  fA = fJacobian_CM * f_Ejectile_Mom_CM_GeV / fPi; // Eqn 21 in paper
 
  // ----------------------------------------------------
  // Jacobian dOmega* / dOmega dimensionless
  // ----------------------------------------------------
  fJacobian_CM_RF  = ( pow((l_Ejectile_rf.Vect()).Mag(),2)*fW) / 
    ( f_Ejectile_Mom_CM * std::abs( ( fProton_Mass + lphoton_rf.E()) * (l_Ejectile_rf.Vect()).Mag() - 
			    ( l_Ejectile_rf.E() * (lphoton_rf.Vect()).Mag() * cos( l_Ejectile_rf.Theta() ) ) ) ); // Differs from next line in photon vect -> lphoton_rf vs r_lphoton
 
  fJacobian_CM_Col = ( ( pow((r_l_Ejectile.Vect()).Mag(),2) * fW ) / // This one is actually used subsequently, so this must be Eqn 20
		       ( f_Ejectile_Mom_CM * std::abs( ( fProton_Mass + r_lphoton.E() ) * (r_l_Ejectile.Vect()).Mag() -
					       ( r_l_Ejectile.E() * (r_lphoton.Vect()).Mag() * cos( r_l_Ejectile.Theta() ) ) ) ) ); 


  //	 cout <<  l_Ejectile_rf.Vect().Mag() << "  " << << << << << << << << endl;
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

  dRootTree->Fill();

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

  fProton_Mom_Col   = fHBeam * 1e3; 
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

  if (rEjectile == "Pi+"){
    total_sig = GetPiPlus_CrossSection(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon);
  }
  else if (rEjectile == "Pi0"){
    total_sig = GetKPlus_CrossSection(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon, rRecoil);
  }
  else if (rEjectile == "K+"){
    total_sig = GetKPlus_CrossSection(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon, rRecoil);
    // SJDK - 31/01/23 - This second case gets the total cross section for the K+ as calculated from Ali's earlier scaling calculation, retain this for comparison.
    //total_sig2 = GetKPlus_CrossSection_Scaling(fT_GeV, fW_GeV, fQsq_GeV, fEpsilon, f_Ejectile_Mass_GeV, rRecoil);
  }
  
  //cout << fT_GeV <<  "  " << fW_GeV << "   " << fQsq_GeV << "   " << "KPlus Paramterisation - " << total_sig << " !!! Old Scaling method - " << total_sig2 << endl;
  
  return total_sig;

}

//// SJDK 21/12/22 - This function needs updating!
//Double_t DEMP_Reaction::GetPi0_CrossSection() {
//
//  double_t sig_total;
//  return sig_total;
//
//}

/*--------------------------------------------------*/
/// Output generator detail
// 06/09/23 SJDK - Formatting of these is all messed up annoyingly, would be nice to get the final numbers to align. They don't currently.
// Cuts are now ordered as they are applied in the generator
void DEMP_Reaction::Detail_Output() {

  DEMPDetails << "Seed used for the Random Number Generator                    " << setw(20) << fSeed         << endl;
  DEMPDetails << endl;
  DEMPDetails << "Total events tried                                           " << setw(20) << fNGenerated   << endl;
  if(UseSolve == true){
    DEMPDetails << "Total events cut                                             " << setw(20) << (qsq_ev + w_ev + w_neg_ev + fNaN + fConserve + t_ev + fNSigmaNeg + fSolveEvents_0Sol) << endl;
    DEMPDetails << "Total events recorded                                        " << setw(20) << fNRecorded    << endl;
    if (fNGenerated != (qsq_ev + w_ev + w_neg_ev + fNaN + fConserve + t_ev + fNSigmaNeg + fNRecorded + fSolveEvents_0Sol)){
      DEMPDetails << "Total events cut + recorded = events tried?                " << setw(20) << "NO! ERROR!" << endl;
    }
    else{
      DEMPDetails << "Total events cut + recorded = events tried?                " << setw(22) << "Yes! :)" << endl;
    }

  }
  else{
    DEMPDetails << "Total events cut                                             " << setw(20) << (qsq_ev + w_ev + w_neg_ev + fNaN + fConserve + t_ev + fNSigmaNeg) << endl;
    DEMPDetails << "Total events recorded                                        " << setw(20) << fNRecorded    << endl;
    if (fNGenerated != (qsq_ev + w_ev + w_neg_ev + fNaN + fConserve + t_ev + fNSigmaNeg + fNRecorded)){
      DEMPDetails << "Total events cut + recorded = events tried?                " << setw(20) << "NO! ERROR!" << endl;
    }
    else{
      DEMPDetails << "Total events cut + recorded = events tried?                " << setw(22) << "Yes! :)" << endl;
    }
  }
  
  DEMPDetails << endl << "Cut details -" << endl;
  DEMPDetails << "Events cut due to qsq < " << fQsq_Min << " or qsq > "<< fQsq_Max << "                        " << setw(20) << qsq_ev << endl;
  DEMPDetails << "Events cut due to negative Wsq value                         " << setw(20) << w_neg_ev      << endl;  
  DEMPDetails << "Events cut due to W < " << fW_Min << " or W > " << fW_Max << "                          " << setw(20) << w_ev << endl;
  if(UseSolve == true){
    DEMPDetails << "Events cut due to solve function finding 0 solutions         " << setw(20) << fSolveEvents_0Sol << endl;
  }
  DEMPDetails << "Events cut due to ejectile (X) energy NaN                    " << setw(20) << fNaN          << endl;
  DEMPDetails << "Events cut due to conservation law check failure             " << setw(20) << fConserve     << endl;
  DEMPDetails << "Events cut due to -t > " << fT_Max << "GeV                      " << setw(30) << t_ev          << endl;
  DEMPDetails << "Events cut due to -ve cross section value                    " << setw(20) << fNSigmaNeg    << endl;

  DEMPDetails << endl << "Conservation law checks details -" << endl;
  DEMPDetails << "Total events PASSING conservation law check with tolerance " << fDiff << setw(17) << conserve   << endl;
  DEMPDetails << "Events cut due to energy conservation check ONLY             " << setw(20) << ene   << endl; 
  DEMPDetails << "Events cut due to momentum conservation check ONLY           " << setw(20) << mom   << endl;
  DEMPDetails << "Events cut due to energy AND momentum conservation checks    " << setw(20) << ene_mom   << endl;
  DEMPDetails << "Events cut due to px conservation law check                  " << setw(20) << mom_px   << endl;
  DEMPDetails << "Events cut due to py conservation law check                  " << setw(20) << mom_py   << endl;
  DEMPDetails << "Events cut due to pz conservation law check                  " << setw(20) << mom_pz   << endl;
  DEMPDetails << "Events cut due to px and py conservation law checks          " << setw(20) << mom_pxpy   << endl;
  DEMPDetails << "Events cut due to px and pz conservation law checks          " << setw(20) << mom_pxpz   << endl;
  DEMPDetails << "Events cut due to py and pz conservation law checks          " << setw(20) << mom_pypz   << endl;
  DEMPDetails << "Events cut due to px, py and pz conservation law checks      " << setw(20) << mom_pxpypz   << endl;

  if(UseSolve == true){
    DEMPDetails << endl << "Solve function, addtional info -" << endl;
    DEMPDetails << "Number of events with 0 Solution                             " << setw(20) << fSolveEvents_0Sol << endl;
    DEMPDetails << "Number of events with 1 Solution                             " << setw(20) << fSolveEvents_1Sol << endl;
    DEMPDetails << "Number of events with 2 Solution                             " << setw(20) << fSolveEvents_2Sol << endl;
  }
  
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
	  << setw(16) << r_l_Ejectile_g.X()
	  << setw(16) << r_l_Ejectile_g.Y()   
	  << setw(16) << r_l_Ejectile_g.Z()  
	  << setw(16) << r_l_Ejectile_g.E()
	  << setw(16) << r_l_Ejectile_g.M() //15/05/23 - Love - Was fX_Mass_GeV
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
	  << setw(16) << r_lscatelecg.M() //15/05/23 - Love- Was fElectron_Mass_GeV
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
	  << setw(16) << l_Recoil_g.X() 
	  << setw(16) << l_Recoil_g.Y()
	  << setw(16) << l_Recoil_g.Z()
	  << setw(16) << l_Recoil_g.E()
	  << setw(16) << l_Recoil_g.M() // 15/05/23 - Love - Was f_Scat_hadron_Mass_GeV
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
	  << setw(14) <<  r_lelectrong.M() // 15/05/23 - Love - Was fElectron_Mass_GeV
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
	 << setw(14) << r_lprotong.M() // 15/05/23 - Love - Was fProton_Mass_GeV
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
	 << setw(14) << r_lscatelecg.M() // 15/05/23 - Love - Was fElectron_Mass_GeV
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

	 << setw(14) << l_Recoil_g.X() 
	 << setw(14) << l_Recoil_g.Y()
	 << setw(14) << l_Recoil_g.Z()
	 << setw(14) << l_Recoil_g.E()
	 << setw(14) <<  l_Recoil_g.M() // 15/05/23 - Love - Was f_Scat_hadron_Mass_GeV
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

	 << setw(14) << r_l_Ejectile_g.X()
	 << setw(14) << r_l_Ejectile_g.Y()   
	 << setw(14) << r_l_Ejectile_g.Z()  
	 << setw(14) << r_l_Ejectile_g.E()
	 << setw(14) <<  r_l_Ejectile_g.M() // 15/05/23 - Love - Was fX_Mass_GeV
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
  
  // HEPMC3 output for Athena/ePIC simulations

  // First line - E - Event# - #Vertices - #Particles
  DEMPOut << std::scientific << std::setprecision(15) << "E" << " "  << print_itt <<  " " << "1" << " " << 5 << endl;
  print_itt++;
  // Second line, Units - U - ENERGY UNIT - DISTANCE UNIT
  DEMPOut << "U" << " " << "GEV" << " " << "MM" << endl;
  // Third line, optional attributes, the weight
  DEMPOut << "A" << " " << "0" << " " << "weight" << " " <<  fEventWeight << endl;
  // Beam particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (4, incoming beam particle)
  DEMPOut << "P" << " " << "1" << " " << "0" << " " << "11" << " " << r_lelectrong.X() << " " << r_lelectrong.Y() << " " << r_lelectrong.Z() << " " << r_lelectrong.E() << " " << r_lelectrong.M() << " " << "4" << endl;
  DEMPOut << "P" << " " << "2" << " " << "0" << " " << "2212" << " " << r_lprotong.X() << " " << r_lprotong.Y() << " " << r_lprotong.Z() << " " << r_lprotong.E() << " " <<  r_lprotong.M()<< " " << "4" << endl;
  // Vertex line - V - 1 - 0 - [1,2]
  DEMPOut << "V" << " " << "-1" << " " << "0" << " " << "[1,2]" << endl;
  // Output particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (1, undecayed physical particle)
  // Scattered electron
  DEMPOut << "P" << " " << "3" << " " << "-1" << " " << "11" << " " << r_lscatelecg.X() << " "  << r_lscatelecg.Y() << " "  << r_lscatelecg.Z() << " " << r_lscatelecg.E() << " " << r_lscatelecg.M() << " " << "1" << endl;
  // Produced meson
  DEMPOut << "P" << " " << "4" << " " << "-1" << " " << PDGtype(produced_X) << " " << r_l_Ejectile_g.X() << " "  << r_l_Ejectile_g.Y() << " "  << r_l_Ejectile_g.Z() << " " << r_l_Ejectile_g.E() << " " << r_l_Ejectile_g.M() << " " << "1" << endl;
  // Recoil hadron
  DEMPOut << "P" << " " << "5" << " " << "-1" << " " << PDGtype(recoil_hadron) << " " << l_Recoil_g.X() << " "  << l_Recoil_g.Y() << " "  << l_Recoil_g.Z() << " " << l_Recoil_g.E() << " " <<  l_Recoil_g.M() << " " << "1" << endl;

}

/*--------------------------------------------------*/ 

bool DEMP_Reaction::SolnCheck()
{
  //
  //  // Double Checking for solution viability
  //  if (TMath::Abs(f_Scat_hadron_Mass-r_l_scat_hadron_solved->M())>1){
  //    //cerr << "Mass Missmatch" << endl;
  //    //cerr << TMath::Abs(proton_mass_mev-Proton->M()) << endl;
  //    return false;
  //  }
  //  if (TMath::Abs(W_in()-W_out())>1){
  //    //cerr << "W Missmatch" << endl;
  //    //cerr << TMath::Abs(W_in()-W_out()) << endl;
  //    return false;
  //  }
  //  *Final = *r_l_scat_hadron_solved + *r_lX_solved;
  //
  //  if (TMath::Abs(Initial->Px()-Final->Px())>1){
  //    //cerr << "Px Missmatch" << endl;
  //    //cerr << TMath::Abs(Initial->Px()-Final->Px()) << endl;
  //    return false;
  //  }
  //
  //  if (TMath::Abs(Initial->Py()-Final->Py())>1){
  //    //cerr << "Py Missmatch" << endl;
  //    //cerr << TMath::Abs(Initial->Py()-Final->Py()) << endl;
  //    return false;
  //  }
  //
  //  if (TMath::Abs(Initial->Pz()-Final->Pz())>1){
  //    //cerr << "Pz Missmatch" << endl;
  //    //cerr << TMath::Abs(Initial->Pz()-Final->Pz()) << endl;
  //    return false;
  //  }
  //
  //  if (TMath::Abs(Initial->E()-Final->E())>1){
  //    return false;
  //  }
  return true;
}

/*--------------------------------------------------*/ 
double DEMP_Reaction::W_in()
{
  return 0;
}

/*--------------------------------------------------*/ 
double DEMP_Reaction::W_out()
{
  return 0;
}

/*--------------------------------------------------*/ 

int DEMP_Reaction::Solve()
{


  VertBeamElec->SetPxPyPzE(r_lelectron.Px(), r_lelectron.Py(), r_lelectron.Pz(), r_lelectron.E());
  VertScatElec->SetPxPyPzE(r_lscatelec.Px(), r_lscatelec.Py(), r_lscatelec.Pz(), r_lscatelec.E());
  Target->SetPxPyPzE(r_lproton.Px(), r_lproton.Py(), r_lproton.Pz(), r_lproton.E());
  *Photon = *VertBeamElec - *VertScatElec;
  *Interaction = *Photon;

  *Initial = *Interaction+*Target;

  theta =  f_Ejectile_Theta_Col;
  phi   =  f_Ejectile_Phi_Col;  

  return this->Solve(theta, phi);
}


int DEMP_Reaction::Solve(double theta, double phi)
{

  W_in_val = W_in();

  if (W_in_val<0){
    return 0;
  }

  UnitVect->SetTheta(theta);
  UnitVect->SetPhi(phi);
  UnitVect->SetMag(1);

  double* pars = new double[9];

  pars[0] = UnitVect->X();
  pars[1] = UnitVect->Y();
  pars[2] = UnitVect->Z();
  pars[3] = Initial->Px();
  pars[4] = Initial->Py();
  pars[5] = Initial->Pz();
  pars[6] = Initial->E();
  pars[7] = f_Ejectile_Mass;
  pars[8] = f_Recoil_Mass;

  F->SetParameters(pars);


  ///*--------------------------------------------------*/ 
  // Looking for the 1st Solution:
  //    If a solution found, then this will be the fist solution. Then we proceed to look for the 2nd solution. 
  //    If no soluion found, then exit solve function

  P = F->GetX(0, 0, pars[6], 0.0001, 10000);

  if (TMath::Abs(F->Eval(P)) < 1){
    fSolveEvents_1Sol++;
  } else {
    fSolveEvents_0Sol++; 
    return 0;
  }

  TLorentzVector * r_l_Ejectile_solved_1_temp = new TLorentzVector();
  TLorentzVector * r_l_Ejectile_solved_2_temp = new TLorentzVector();

  Float_t r_l_Ejectile_E = sqrt( pow(P*pars[0],2) + pow(P*pars[1],2) + pow(P*pars[2],2) + pow(f_Ejectile_Mass,2) );
  r_l_Ejectile_solved_1_temp->SetPxPyPzE(P*pars[0], P*pars[1],  P*pars[2], r_l_Ejectile_E);

  ///*--------------------------------------------------*/ 
  // Looking for the 2nd Solution 

  P2 = F->GetX(0, P+100, pars[6], 0.0001, 10000);
  Float_t r_l_Ejectile_E_2 = sqrt( pow(P2 * pars[0],2) + pow(P2 * pars[1],2) + pow(P2 * pars[2],2) + pow(f_Ejectile_Mass,2) );
  r_l_Ejectile_solved_2_temp->SetPxPyPzE(P2 * pars[0], P2 * pars[1],  P2 * pars[2], r_l_Ejectile_E_2);

  ///*--------------------------------------------------*/ 
  // If a valid 2nd solution is found, then we are certian that there are two solutions.
  //   - We then increament the counter for 2nd solution scenario
  //   - We then decreament the counter for the 1st solution scenario 

  if (TMath::Abs(F->Eval(P2)) < 1){
    fSolveEvents_2Sol++;
    fSolveEvents_1Sol--;
    if ( Int_t(CoinToss->Uniform(0,100)) < 50) {
      r_l_Ejectile_solved.SetPxPyPzE(r_l_Ejectile_solved_1_temp->X(), r_l_Ejectile_solved_1_temp->Y(), r_l_Ejectile_solved_1_temp->Z(), r_l_Ejectile_solved_1_temp->E());
    } else {
      r_l_Ejectile_solved.SetPxPyPzE(r_l_Ejectile_solved_2_temp->X(), r_l_Ejectile_solved_2_temp->Y(), r_l_Ejectile_solved_2_temp->Z(), r_l_Ejectile_solved_2_temp->E());
    }
  }
  else {
    r_l_Ejectile_solved.SetPxPyPzE(r_l_Ejectile_solved_1_temp->X(), r_l_Ejectile_solved_1_temp->Y(), r_l_Ejectile_solved_1_temp->Z(), r_l_Ejectile_solved_1_temp->E());
  }

  ///*--------------------------------------------------*/ 
  /// Solve for the recoil information with the "solved" Ejectile informaiton
  TLorentzVector * r_l_hadron_temp= new TLorentzVector();
  *r_l_hadron_temp = *Initial- r_l_Ejectile_solved;
  r_l_Recoil_solved.SetPxPyPzE(r_l_hadron_temp->Px(), r_l_hadron_temp->Py(), r_l_hadron_temp->Pz(), r_l_hadron_temp->E());

  delete r_l_Ejectile_solved_1_temp;
  delete r_l_Ejectile_solved_2_temp;

  delete r_l_hadron_temp;
  delete[] pars;
 
  return 1;

}
