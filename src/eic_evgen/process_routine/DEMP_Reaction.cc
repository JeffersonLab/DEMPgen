#include "reaction_routine.h"
#include "eic.h"

#include <sys/stat.h>

using namespace std;

DEMP_Reaction::DEMP_Reaction(){ 

  cout << "Program Start" << endl;

}

//-------------------Love Preet - Added for phase space factor calculations----------------------------//
 
Double_t DEMP_Reaction::calculate_psf_t(const TLorentzVector& psf_photon, const TLorentzVector& psf_ejectile){ // Calculate the Lorentz variable (t) for Phase space factor (psf)/
  return -1.0*((psf_photon - psf_ejectile).Mag2());
}
  
void DEMP_Reaction::calculate_psf_max_min(double value, double& maxValue, double& minValue){// Find the max and min values for scattered electron's energy, theta and ejetile's theta
  if (value > maxValue) {
    maxValue = value;
  }
  if (value < minValue) {
    minValue = value;
  }
}

//-------------------Love Preet - Added for phase space factor (PSF) calculations----------------------------//
Double_t DEMP_Reaction::psf(){

  // Added a print out to clarify what is happening
  cout << "Beginning phase space calculation -" << endl;

  //----Note that all the calculations have been done in GeV for the psf...............................................//
  //-----------------------------Calculate variables for scattered electron...........................................//
  psf_ScatElec_E_Stepsize = (fEBeam*(fScatElec_E_Hi - fScatElec_E_Lo))/psf_steps; //in GeV // Defined stepsize for scattered electron's energy, theta, and phi 
  psf_ScatElec_Theta_Stepsize = (fScatElec_Theta_F - fScatElec_Theta_I)/psf_steps; // in radian
  psf_ScatElec_Phi_Stepsize = (2.0 * fPi)/4.0;
   
  psf_Ejec_Theta_Stepsize = (fEjectileX_Theta_F - fEjectileX_Theta_I)/psf_steps; // in radian  // Defined stepsize for ejectile theta
   
  for (int psf_E = 0; psf_E <= psf_steps; psf_E++){ // Loop over the scattered electron's energy
    // SJDK - 29/04/24 - Added progress report
    dFractTime = time(0); 

    if ( psf_E % ( psf_steps / 10 ) == 0 ){
      cout << ((1.0*psf_E)/(1.0*psf_steps))*100.0 << setw(4) << " % of phase space checked"  
	   << "   Day: " <<  dFractTime.GetDay() 
	   << "   Time:   " << dFractTime.GetHour() 
	   << ":" << dFractTime.GetMinute() 
	   << ":" << dFractTime.GetSecond() 
	   << endl;	  
    }

    psf_ScatElec_E = (fEBeam * fScatElec_E_Lo + (psf_E * psf_ScatElec_E_Stepsize));  
    psf_ScalElec_Mom = sqrt(pow(psf_ScatElec_E,2) - pow(fElectron_Mass_GeV,2));

    for (int psf_Theta = 0; psf_Theta <= psf_steps; psf_Theta++){  // Loop over the scattered electron's theta

      psf_ScatElec_Theta = (fScatElec_Theta_I + (psf_Theta * psf_ScatElec_Theta_Stepsize));

      for (int psf_Phi = 0; psf_Phi <= 4; psf_Phi++){ // Loop over the scattered electron's Phi

	psf_ScatElec_Phi = (0.0 + (psf_Phi * psf_ScatElec_Phi_Stepsize));
	psf_scatelec.SetPxPyPzE(psf_ScalElec_Mom * sin(psf_ScatElec_Theta) * cos( psf_ScatElec_Phi),   // Scattered's electron four momentum
				psf_ScalElec_Mom * sin(psf_ScatElec_Theta) * sin( psf_ScatElec_Phi),
				psf_ScalElec_Mom * cos(psf_ScatElec_Theta),
				psf_ScatElec_E);
	psf_photon = r_lelectrong - psf_scatelec; // virtual photon four momentum
	psf_Q2 = -1.*(psf_photon.Mag2());         // Lorentz variable Q2
	psf_W  =  ((psf_photon + r_lprotong).Mag()); // Lorentz variable W
	psf_W2 =  ((psf_photon + r_lprotong).Mag2()); // Lorentz variable W2
	
	if (( psf_Q2 >= fQsq_Min && psf_Q2 <= fQsq_Max) && psf_W2 >= 0.0 && ( psf_W >= fW_Min && psf_W <= fW_Max)){
	  for (int psf_Ejec_Theta = 0; psf_Ejec_Theta <= psf_steps; psf_Ejec_Theta++){ // Loop over the scattered ejectile's theta
  
	    //-----------------------------Calculate variables for ejectile..........................................//      
	    psf_Ejectile_Theta = (fEjectileX_Theta_I + (psf_Ejec_Theta *psf_Ejec_Theta_Stepsize));
	    psf_Ejectile_Phi = 0.0;
 
	    double fupx = sin( psf_Ejectile_Theta ) * cos( psf_Ejectile_Phi );
	    double fupy = sin( psf_Ejectile_Theta ) * sin( psf_Ejectile_Phi );
	    double fupz = cos( psf_Ejectile_Theta );
  
	    double fuqx = sin( psf_photon.Theta() ) * cos( psf_photon.Phi() );
	    double fuqy = sin( psf_photon.Theta() ) * sin( psf_photon.Phi() );
	    double fuqz = cos( psf_photon.Theta() );

	    double fa = -(psf_photon.Vect()).Mag() * ( fupx * fuqx +  fupy * fuqy +  fupz * fuqz);
	    double fb = pow ( (psf_photon.Vect()).Mag() , 2 );
	    double fc = psf_photon.E() +  fProton_Mass_GeV;

	    fa = ( fa - std::abs( (r_lprotong.Vect()).Mag() ) * ( ( ( r_lprotong.X() / (r_lprotong.Vect()).Mag() ) * fupx ) + 
								  ( ( r_lprotong.Y() / (r_lprotong.Vect()).Mag() ) * fupy ) + 
								  ( ( r_lprotong.Z() / (r_lprotong.Vect()).Mag() ) * fupz ) ) );
   
	    double factor = ( pow( (r_lprotong.Vect()).Mag() , 2 ) + 2.0 * (psf_photon.Vect()).Mag() * (r_lprotong.Vect()).Mag() *
			      ( ( ( r_lprotong.X() / (r_lprotong.Vect()).Mag() ) * fuqx ) + 
				( ( r_lprotong.Y() / (r_lprotong.Vect()).Mag() ) * fuqy ) + 
				( ( r_lprotong.Z() / (r_lprotong.Vect()).Mag() ) * fuqz ) ) );
  
	    fb =  fb + factor;
	    fc = psf_photon.E() + r_lprotong.E();

	    double ft = fc * fc - fb + f_Ejectile_Mass_GeV * f_Ejectile_Mass_GeV -  f_Recoil_Mass_GeV *  f_Recoil_Mass_GeV;
     
	    double fQA = 4.0 * ( fa * fa - fc * fc );
	    double fQB = 4.0 * fc * ft;

	    double fQC = -4.0 * fa * fa * f_Ejectile_Mass_GeV * f_Ejectile_Mass_GeV - ft * ft;
 
	    fradical = fQB * fQB - 4.0 * fQA * fQC;
	    
	    fepi1 = ( -fQB - sqrt( fradical ) ) / ( 2.0 * fQA );
	    fepi2 = ( -fQB + sqrt( fradical ) ) / ( 2.0 * fQA );
 
	    psf_ejectile.SetPxPyPzE( ( sqrt( pow( fepi1 , 2) - pow(f_Ejectile_Mass_GeV , 2) ) ) * sin( psf_Ejectile_Theta ) * cos(  psf_Ejectile_Phi ),
				     ( sqrt( pow( fepi1 , 2) - pow(f_Ejectile_Mass_GeV , 2) ) ) * sin( psf_Ejectile_Theta ) * sin(  psf_Ejectile_Phi ),
				     ( sqrt( pow( fepi1 , 2) - pow(f_Ejectile_Mass_GeV , 2) ) ) * cos( psf_Ejectile_Theta ),
				     fepi1 );

	    if (!TMath::IsNaN(psf_ejectile.E())) { // If the ejectile's energy is a number
   
	      psf_t = calculate_psf_t(psf_photon, psf_ejectile);
  
	      if (psf_t >= 0.0 && psf_t <= fT_Max) { // If the t value is in this range
 
		calculate_psf_max_min(psf_ScatElec_E, psf_ScatElec_E_max, psf_ScatElec_E_min);
		calculate_psf_max_min(psf_ScatElec_Theta, psf_ScatElec_Theta_max, psf_ScatElec_Theta_min);
		calculate_psf_max_min(psf_Ejectile_Theta, psf_Ejectile_Theta_max, psf_Ejectile_Theta_min);
  
	      } // if statement over psf_t
	    } // if statement over ejectile energy
	  } // End of for loop over psf_Theta_Ejec
	} // If condition over psf_Q2,psf_W,psf_W2
      } // End of for loop over psf_Phi
    } // End of for loop over psf_Theta
  } // End of for loop over psf_Phi
    
  if ((psf_ScatElec_E_max == psf_ScatElec_Theta_max) && (psf_ScatElec_Theta_max == psf_Ejectile_Theta_max) && (psf_Ejectile_Theta_max == psf_ScatElec_E_max)){
    cout << "!!!!! - ERROR - !!!!!" << endl;
    cout << "Specified scattered electron and ejectile energy/angle ranges in input .json file did not yield any valid allowed phase space." << endl << "Processing stopped, please re-specify values in input .json file." << endl;
    cout << "!!!!! - ERROR - !!!!!" << endl;
    fPSF = 0.0;
  }
  else if ((psf_ScatElec_E_min == psf_ScatElec_Theta_min) && (psf_ScatElec_Theta_min == psf_Ejectile_Theta_min) && (psf_Ejectile_Theta_min == psf_ScatElec_E_min)){
    cout << "!!!!! - ERROR - !!!!!" << endl;
    cout << "Specified scattered electron and ejectile energy/angle ranges in input .json file did not yield any valid allowed phase space." << endl << " Processing stopped, please re-specify values in input .json file." << endl;
    cout << "!!!!! - ERROR - !!!!!" << endl;
    fPSF = 0.0;
  } 
  else if ((psf_ScatElec_E_max == psf_ScatElec_E_min) || (psf_ScatElec_Theta_max == psf_ScatElec_Theta_min) || (psf_Ejectile_Theta_max == psf_Ejectile_Theta_min)) {
    cout << "!!!!! - ERROR - !!!!!" << endl;
    cout << "Specified scattered electron and ejectile energy/angle ranges in input .json file include one (or more) values with equal min/max values in range." << endl << "Processing stopped, please re-specify values in input .json file." << endl;
    cout << "!!!!! - ERROR - !!!!!" << endl;
    fPSF = 0.0;
  } 
  else {
    // Calculate the PSF by expanding the limits on both sides by a step size of 1.0.
    psf_ScatElec_E_min      = psf_ScatElec_E_min     - ( 1.0 * psf_ScatElec_E_Stepsize ); // Units are in GeV
    psf_ScatElec_E_max      = psf_ScatElec_E_max     + ( 1.0 * psf_ScatElec_E_Stepsize );
    psf_ScatElec_Theta_min  = psf_ScatElec_Theta_min - ( 1.0 * psf_ScatElec_Theta_Stepsize ); // Units are in Radian
    psf_ScatElec_Theta_max  = psf_ScatElec_Theta_max + ( 1.0 * psf_ScatElec_Theta_Stepsize );
    psf_Ejectile_Theta_min  = psf_Ejectile_Theta_min - ( 1.0 * psf_Ejec_Theta_Stepsize ); // Units are in Radian
    psf_Ejectile_Theta_max  = psf_Ejectile_Theta_max + ( 1.0 * psf_Ejec_Theta_Stepsize );
  	
    if ( psf_ScatElec_E_min < ( fEBeam * fScatElec_E_Lo ) || psf_ScatElec_E_min == ( fEBeam * fScatElec_E_Lo ) ) { psf_ScatElec_E_min = ( fEBeam * fScatElec_E_Lo ); }
    if ( psf_ScatElec_E_max > ( fEBeam * fScatElec_E_Hi ) || psf_ScatElec_E_max == ( fEBeam * fScatElec_E_Hi ) ) { psf_ScatElec_E_max = ( fEBeam * fScatElec_E_Hi ); }
    if ( psf_ScatElec_Theta_min < fScatElec_Theta_I || psf_ScatElec_Theta_min == fScatElec_Theta_I ) { psf_ScatElec_Theta_min = fScatElec_Theta_I; }
    if ( psf_ScatElec_Theta_max > fScatElec_Theta_F || psf_ScatElec_Theta_max == fScatElec_Theta_F ) { psf_ScatElec_Theta_max = fScatElec_Theta_F; }
    if ( psf_Ejectile_Theta_min < fEjectileX_Theta_I || psf_Ejectile_Theta_min == fEjectileX_Theta_I ) { psf_Ejectile_Theta_min = fEjectileX_Theta_I; }
    if ( psf_Ejectile_Theta_max > fEjectileX_Theta_F || psf_Ejectile_Theta_max == fEjectileX_Theta_F ) { psf_Ejectile_Theta_max = fEjectileX_Theta_F; }
  
    fPSF = (( psf_ScatElec_E_max - psf_ScatElec_E_min ) *( cos( psf_ScatElec_Theta_max ) - cos( psf_ScatElec_Theta_min ) ) * 2 * fPI *( cos( psf_Ejectile_Theta_max ) - cos( psf_Ejectile_Theta_min ) ) * 2 * fPI );
  }
  return fPSF;
}

/*--------------------------------------------------*/
/// DEMP_Reaction
DEMP_Reaction::DEMP_Reaction(TString particle_str, TString hadron_str){ 

  rEjectile = particle_str;
  rRecoil = hadron_str;

}

DEMP_Reaction::~DEMP_Reaction(){

  // Text output data file and generation information 
  DEMPOut.close();
  DEMPDetails.close();

  if (gROOTOut == true){
    // Diagnostic root file with plots
    dRootFile->Write(); // Write the contents of the ROOT file to disk.
    delete dRootTree;   // Delete the dynamically allocated memory for the ROOT tree.
    dRootFile->Close(); // Close the ROOT file.
    delete dRootFile;   // Delete the dynamically allocated memory for the ROOT file object.
  }
}

void DEMP_Reaction::process_reaction(){
 
  Init();
  
  if (fPSF <= 0){// If phase space factor is zero or less than zero, stop further processing // Love Preet - Added for phase space factor calculations
    return;
  }
  
  //-------------------Love Preet - Added to print out the considered ranges for PSF calculations----------------------------//
  
  cout << "------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "Specified scattered electron and ejectile energy/angle ranges in input .json file are ajusted to ensure that they fall within the allowed phase space." << endl;
  cout << "Ee_Low = " << psf_ScatElec_E_min << ", Ee_High = " << psf_ScatElec_E_max << ", e_Theta_Low = " << psf_ScatElec_Theta_min * TMath::RadToDeg() << ", e_Theta_High = " << psf_ScatElec_Theta_max * TMath::RadToDeg() << ", EjectileX_Theta_Low = " << psf_Ejectile_Theta_min * TMath::RadToDeg() << ", EjectileX_Theta_High = " << psf_Ejectile_Theta_max * TMath::RadToDeg() << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  
  //-------------------Love Preet - Added to print out the considered ranges for PSF calculations----------------------------//
  
  if (gOutputType == "Pythia6"){
    DEMPReact_Pythia6_Out_Init();
  }
  else if (gOutputType == "HEPMC3"){
    DEMPReact_HEPMC3_Out_Init();
  }

  for( long long int i = 0; i < rNEvents; i++ ){
 
    rNEvent_itt = i;
    fNGenerated ++;
 
    Progress_Report();  // This happens at each 10% of the total event is processed
    Processing_Event();
  }

  //-------------------Love Preet - Added for actual phase space factor calculations----------------------------//
  fPSF_org = (( fScatElec_Energy_Col_max * fm - fScatElec_Energy_Col_min * fm ) *( cos( fScatElec_Theta_Col_max ) - cos( fScatElec_Theta_Col_min ) ) * 2 * fPI *( cos( f_Ejectile_Theta_Col_max ) - cos( f_Ejectile_Theta_Col_min ) ) * 2 * fPI ); // Calculate the actual phase space factor, this is determined from the actual scattered electron/ejectile values in the generated events that passed all cuts
  
  //-------------------Love Preet - Added for actual phase space factor calculations----------------------------// 

  Detail_Output();

}

void DEMP_Reaction::Init(){

  pim* myPim;

  pd = dynamic_cast<pim*>(myPim);
	
  rEjectile_charge = ExtractCharge(rEjectile);

  const char* dir_name = "OutputFiles";
  struct stat sb;

  if (stat(dir_name, &sb) == 0) {
    cout << "Output file directory found from DEMPgen directory - " << dir_name  << endl;
  }
  else {
    cout << "Output file directory not found from DEMPgen directory - " << dir_name  << endl;
    cout << "Making OutputFiles directory!" << endl;
    mkdir(dir_name,0777);
  } 
  
  sTFile = Form("./%s/eic_%s.txt", dir_name, gfile_name.Data());
  sLFile = Form("./%s/eic_input_%s.dat", dir_name, gfile_name.Data());

  DEMPOut.open( sLFile.c_str() );
  DEMPDetails.open( sTFile.c_str() );

  if (gROOTOut == true){ // Only initialise and open root file if output is enabled
    sDFile = Form("./%s/eic_%s.root", dir_name, gfile_name.Data()); // LovePreet changed to make the files name consistent
    dRootFile = new TFile(sDFile.c_str(),"RECREATE"); 
    dRootTree = new TTree("Events", "Description of a tree");  //Love Preet added all these new braches to be stored in a root ttree
    dRootTree->Branch("EventWeight", &fEventWeight, "EventWeight/D");
    dRootTree->Branch("scat_e_px", &scat_e_px, "scat_e_px/D");
    dRootTree->Branch("scat_e_py", &scat_e_py, "scat_e_py/D");
    dRootTree->Branch("scat_e_pz", &scat_e_pz, "scat_e_pz/D");
    dRootTree->Branch("scat_e_E",  &scat_e_E, "scat_e_E/D");
    dRootTree->Branch("ejec_px", &ejec_px, "ejec_px/D");
    dRootTree->Branch("ejec_py", &ejec_py, "ejec_py/D");
    dRootTree->Branch("ejec_pz", &ejec_pz, "ejec_pz/D");
    dRootTree->Branch("ejec_E",  &ejec_E, "ejec_E/D");
    dRootTree->Branch("rclH_px", &rclH_px, "rclH_px/D");
    dRootTree->Branch("rclH_py", &rclH_py, "rclH_py/D");
    dRootTree->Branch("rclH_pz", &rclH_pz, "rclH_pz/D");
    dRootTree->Branch("rclH_E",  &rclH_E, "rclH_E/D"); 
    dRootTree->Branch("Q2", &fQsq_GeV, "Q2/D");
    dRootTree->Branch("W",  &fW_GeV, "W/D");
    dRootTree->Branch("t", &fT_GeV, "t/D");
    dRootTree->Branch("x_b", &fx, "x_b/D");
    dRootTree->Branch("y_E", &fy, "y_E/D");
  }
  /*--------------------------------------------------*/
  qsq_ev = 0, t_ev = 0, w_neg_ev = 0, w_ev = 0;
  rNEvents = fNEvents;
  rNEvent_itt = 0;
  
  // 02/06/21 SJDK 
  // Set these values once the beam energies are read in
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
  cout << "Define: " << fElectron_MomZ_Col << "    "<< fElectron_Mom_Col << "  " << cos(fElectron_Theta_Col) << endl;

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
  
  fPSF = psf(); // Love Preet - Added for phase space factor calculations
  if (fPSF <= 0){ // if phase space factor is zero or less than zero, stop further processing // Love Preet - Added for phase space factor calculations
    return;
  }
   
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

void DEMP_Reaction::Processing_Event(){

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
  //fScatElec_Theta_Col  = acos( fRandom->Uniform( cos( fScatElec_Theta_I ) , cos( fScatElec_Theta_F ) ) );
  fScatElec_Theta_Col  = acos( fRandom->Uniform( cos( psf_ScatElec_Theta_min ) , cos( psf_ScatElec_Theta_max ) ) ); // Love Preet changed to the allowed region calculated from psf 
  fScatElec_Phi_Col    = fRandom->Uniform( 0 , 2.0 * fPi);
  //fScatElec_Energy_Col = fRandom->Uniform( fScatElec_E_Lo * fElectron_Energy_Col , fScatElec_E_Hi * fElectron_Energy_Col ); // in MeV
  fScatElec_Energy_Col = fRandom->Uniform( psf_ScatElec_E_min * fK, psf_ScatElec_E_max * fK ); // in MeV // Love Preet changed to the allowed region calculated from psf 
 
  // ----------------------------------------------------
  // Produced ejectile in Collider frame
  // ----------------------------------------------------  
  //f_Ejectile_Theta_Col      = acos( fRandom->Uniform( cos(f_Ejectile_Theta_I), cos(f_Ejectile_Theta_F ) ) ); 
  f_Ejectile_Theta_Col      = acos( fRandom->Uniform( cos( psf_Ejectile_Theta_min ), cos( psf_Ejectile_Theta_max ) ) ); // Love Preet changed to the allowed region calculated from psf 
  f_Ejectile_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi );  
  
  //---------------------------------------------------------------------
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
  
  // ----------------------------------------------------
  // W square, Invariant Mass (P_g + P_p)^2
  // ----------------------------------------------------
 
  TLorentzVector lwg;
  lwg = r_lprotong + r_lphotong;
  fW_GeV    = lwg.Mag();
  fWSq_GeV  = lwg.Mag2();
        
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
  //*--------------------------------------------------*/ 
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
  
  if( pd->CheckLaws(r_lelectron, r_lproton, r_lscatelec, r_l_Ejectile, l_Recoil) !=1 ){ // Love Preet changed the order of the conditions so that can identify unphysical events based on NaN ejectile energy, conservation law, and negative W.
    fConserve++;
    return;
  }	
   
  if ( fWSq_GeV < 0 ) { 
    w_neg_ev++;
    return;
  } 
  
  // SJDK 03/04/23 - Qsq an W ranges now variables set by particle type in .json read in. See eic.cc
  if ( fQsq_GeV < fQsq_Min || fQsq_GeV > fQsq_Max ) {
    qsq_ev++;
    return;
  }
  
  if ( fW_GeV < fW_Min || fW_GeV > fW_Max ) { // SJDK 03/04/23 - Switched to the new variable, set by particle type
    w_ev++;
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

  //fBeta_CM_RF        = (lphoton_rf.Vect()).Mag() / (lphoton_rf.E() + r_lhadron_beam_mass);
  fBeta_CM_RF        = ( lphoton_rfg.Vect() ).Mag() / ( lphoton_rfg.E() + ( r_lhadron_beam_mass * fm ) ); // involved quantities are in GeV

  //fGamma_CM_RF       = (lphoton_rf.E() + r_lhadron_beam_mass) / fW;
  fGamma_CM_RF       = ( lphoton_rfg.E() + ( r_lhadron_beam_mass * fm ) ) / ( fW * fm ); // involved quantities are in GeV
  f_Ejectile_Energy_CM       = (pow(fW, 2) + pow(f_Ejectile_Mass,2) - pow(f_Recoil_Mass,2) ) / (2.0* fW);    
  f_Ejectile_Mom_CM          = sqrt(pow(f_Ejectile_Energy_CM,2) - pow(f_Ejectile_Mass,2));    
  //f_Ejectile_Energy_CM_GeV   = f_Ejectile_Energy_CM / 1000.0;
  //f_Ejectile_Mom_CM_GeV      = f_Ejectile_Mom_CM / 1000.0;
  f_Ejectile_Energy_CM_GeV   = f_Ejectile_Energy_CM * fm;
  f_Ejectile_Mom_CM_GeV      = f_Ejectile_Mom_CM * fm;

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
 
  //double fTheta_EEp = (lelectron_rf.Vect()).Angle(lscatelec_rf.Vect());
  double fTheta_EEp = ( lelectron_rfg.Vect() ).Angle( lscatelec_rfg.Vect() ); // involved quantities are in GeV

  fEpsilon = 1.0 / ( 1.0 + 2.0 * ( pow( (lphoton_rfg.Vect()).Mag(),2)/fQsq_GeV ) * pow( tan( fTheta_EEp / 2 ) , 2 ) );

  // ----------------------------------------------------
  // Virtual Photon flux factor in units of 1/(GeV*Sr)
  // ----------------------------------------------------
  // fFlux_Factor_Col = (fAlpha/(2.0*pow(fPi,2))) * (r_lscatelecg.E() / r_lelectrong.E()) * 
  //   ( pow(fW_GeV,2) - pow(fProton_Mass_GeV,2) ) / (2.0*fProton_Mass_GeV*fQsq_GeV*(1.0 - fEpsilon));
  
  fFlux_Factor_Col = ( fAlpha / ( 2.0 * pow( fPi , 2 ) ) ) * ( r_lscatelecg.E() / r_lelectrong.E() ) *  // All are in GeV (more organised form)
    ( pow( fW_GeV , 2 ) - pow( fProton_Mass_GeV , 2 ) ) / 
    ( 2.0 * fProton_Mass_GeV * fQsq_GeV * ( 1.0 - fEpsilon ) );
         
  fFlux_Factor_RF = ( fAlpha / ( 2.0 * pow( fPi , 2 ) ) ) * ( lscatelec_rfg.E() / lelectron_rfg.E() ) *
    ( pow( fW_GeV , 2 ) - pow( fProton_Mass_GeV , 2 ) ) /
    ( 2.0 * fProton_Mass_GeV * fQsq_GeV * ( 1.0 - fEpsilon ) );
    
  // cout<<"  "<<fFlux_Factor_Col<<"  "<<fFlux_Factor_RF<<endl;

  // ----------------------------------------------------
  //  Jacobian  dt/dcos(theta*)dphi in units of GeV2/sr
  // ----------------------------------------------------
  // fJacobian_CM = ( (lphoton_rfg.Vect()).Mag() - fBeta_CM_RF * lphoton_rfg.E() ) / ( fGamma_CM_RF * ( 1.0 - pow(fBeta_CM_RF,2) ) ); // Eqn 22 in paper
  fJacobian_CM = fGamma_CM_RF * ( ( lphoton_rfg.Vect() ).Mag() - ( fBeta_CM_RF * lphoton_rfg.E() ) ); // All are in GeV // Love Preeet changed it 
  
  fA = fJacobian_CM * f_Ejectile_Mom_CM_GeV / fPi; // Eqn 21 in paper
 
  // ----------------------------------------------------
  // Jacobian dOmega* / dOmega dimensionless
  // ----------------------------------------------------
  /*fJacobian_CM_RF  = ( pow((l_Ejectile_rf.Vect()).Mag(),2)*fW) / 
    ( f_Ejectile_Mom_CM * std::abs( ( fProton_Mass + lphoton_rf.E()) * (l_Ejectile_rf.Vect()).Mag() - 
    ( l_Ejectile_rf.E() * (lphoton_rf.Vect()).Mag() * cos( l_Ejectile_rf.Theta() ) ) ) ); // Differs from next line in photon vect -> lphoton_rf vs r_lphoton
 
    fJacobian_CM_Col = ( ( pow((r_l_Ejectile.Vect()).Mag(),2) * fW ) / // This one is actually used subsequently, so this must be Eqn 20
    ( f_Ejectile_Mom_CM * std::abs( ( fProton_Mass + r_lphoton.E() ) * (r_l_Ejectile.Vect()).Mag() -
    ( r_l_Ejectile.E() * (r_lphoton.Vect()).Mag() * cos( r_l_Ejectile.Theta() ) ) ) ) ); */

  fBeta_Col  = ( r_lphotong.Vect() + r_lprotong.Vect() ).Mag() / ( r_lphotong.E() + r_lprotong.E() ); // All are in GeV // Love Preeet changed it 
  fGamma_Col = ( r_lphotong.E() + r_lprotong.E() ) / ( fW * fm );
  ftheta_Col = ( r_lphotong.Angle( r_l_Ejectile_g.Vect() ) ); // angle between the virtualphoton and the ejectile in the collider frame
   
  fJacobian_CM_Col = ( pow( ( r_l_Ejectile_g.Vect()).Mag() , 2 ) /
		       fGamma_Col * f_Ejectile_Mom_CM_GeV * ( ( r_l_Ejectile_g.Vect() ).Mag() - ( fBeta_Col * r_l_Ejectile_g.E() * cos( ftheta_Col ) ) ) );


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
  // fSigma_Col = r_fSig * fFlux_Factor_RF * fA * fJacobian_CM_Col; // Love Preet changed flux factor from collider to rest frame

  if ( ( fSigma_Col <= 0 ) || std::isnan( fSigma_Col ) ) { 
    fNSigmaNeg ++;
    return;
  }
     
  // -----------------------------------------------------------------------------------------------------------
  //             Lab cross section     Phase Space   Conversion     Luminosity                Total events tried
  // Hz        = ub / ( sr^2 * GeV ) * GeV * sr^2 * ( cm^2 / ub ) * ( # / ( cm^2 * sec ) ) / ( # )

  // SJDK 24/06/21 - Explicitly taking the absolute value of the weight such that the value is positive! Shouldn't matter since any -ve cross section events should be dumped above
  //fEventWeight = abs(fSigma_Col * fPSF * fuBcm2 * fLumi / fNEvents);   // in Hz
  fEventWeight = fSigma_Col * fPSF * fuBcm2 * fLumi / fNEvents;   // in Hz // Love Preet removed the abs on the fEventWeight
  
  if ( ( fEventWeight <= 0 ) || std::isnan( fEventWeight ) ) {  // Love Preet added new counter to track negative and NaN weights
    fNWeightNeg ++;
    return;
  }
  
  fNRecorded++;
  fRatio = fNRecorded / fNGenerated;
  //Love Preet - Added for actual phase space factor calculations 
  calculate_psf_max_min( fScatElec_Energy_Col, fScatElec_Energy_Col_max, fScatElec_Energy_Col_min ); // -> Love Preet added to find the max and min values to calculate the actual PSF
  calculate_psf_max_min( fScatElec_Theta_Col,  fScatElec_Theta_Col_max,  fScatElec_Theta_Col_min );
  calculate_psf_max_min( f_Ejectile_Theta_Col, f_Ejectile_Theta_Col_max, f_Ejectile_Theta_Col_min ); 
   
  if (gOutputType == "Pythia6"){
    DEMPReact_Pythia6_Output();
  }
  else if (gOutputType == "LUND"){
    Lund_Output();
  }
  else if (gOutputType == "HEPMC3"){
    DEMPReact_HEPMC3_Output();
  }

  if (gROOTOut == true){
    scat_e_px =  r_lscatelecg.X(); // Love Preet - Added to be stored in the root tree
    scat_e_py =  r_lscatelecg.Y();
    scat_e_pz =  r_lscatelecg.Z();
    scat_e_E  =  r_lscatelecg.E();
    ejec_px   =  r_l_Ejectile_g.X();
    ejec_py   =  r_l_Ejectile_g.Y();
    ejec_pz   =  r_l_Ejectile_g.Z();
    ejec_E    =  r_l_Ejectile_g.E(); 
    rclH_px   =  l_Recoil_g.X();
    rclH_py   =  l_Recoil_g.Y();
    rclH_pz   =  l_Recoil_g.Z();
    rclH_E    =  l_Recoil_g.E();
    dRootTree->Fill();
  }
}

void DEMP_Reaction::Progress_Report(){

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

TLorentzVector DEMP_Reaction::GetProtonVector_lab(){

  // Crossing angle
  //	 fProton_Theta_Col = 0.050;
  //	 fProton_Theta_Col = 0.025;
  // SJDK - 12/01/22
  // Set crossing angle to 0 for fun4all, also required for ATHENA simulations
  fProton_Theta_Col = 0.0;

  ///*--------------------------------------------------*/
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

void DEMP_Reaction::Consider_Proton_Fermi_Momentum(){

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

TLorentzVector DEMP_Reaction::GetElectronVector_lab(){

  fElectron_Energy_Col = fElectron_Kin_Col; 
  fElectron_Mom_Col    = sqrt( pow(fElectron_Energy_Col , 2) - pow(fElectron_Mass , 2) );
  fElectron_Theta_Col  = fPi;
  fElectron_Phi_Col    = 0.0;
  fElectron_MomZ_Col   = fElectron_Mom_Col * cos(fElectron_Theta_Col);  
  fElectron_MomX_Col   = fElectron_Mom_Col * sin(fElectron_Theta_Col) * cos(fElectron_Phi_Col);
  fElectron_MomY_Col   = fElectron_Mom_Col * sin(fElectron_Theta_Col) * sin(fElectron_Phi_Col);  

  //cout << "Define: " << fElectron_MomZ_Col << "    "<< fElectron_Mom_Col << "  " << cos(fElectron_Theta_Col) << endl;  
        
  TLorentzVector  lelectron( fElectron_MomX_Col, fElectron_MomY_Col, fElectron_MomZ_Col, fElectron_Energy_Col);

  return lelectron;	 

}

Double_t DEMP_Reaction::Get_Phi_X_LeptonPlane_RF(){

  fCos_Phi_X_LeptonPlane_RF = ( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ); // hep-ph/0410050v2
  fSin_Phi_X_LeptonPlane_RF = ( ( v3LxP.Dot( v3PhotonUnit  ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ); // hep-ph/0410050v2    
  if ( fSin_Phi_X_LeptonPlane_RF >= 0 )
    fPhi_X_LeptonPlane_RF   = fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) );
  if ( fSin_Phi_X_LeptonPlane_RF < 0 )
    fPhi_X_LeptonPlane_RF   = 360.0 - std::abs( fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxP ) ) / ( v3QUnitxL.Mag() * v3QUnitxP.Mag() ) ) );

  return fPhi_X_LeptonPlane_RF;

}

Double_t DEMP_Reaction::Get_Phi_TargPol_LeptonPlane_RF(){

  fCos_Phi_TargPol_LeptonPlane_RF = ( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ); // hep-ph/0410050v2
  fSin_Phi_TargPol_LeptonPlane_RF = ( ( v3LxS.Dot( v3PhotonUnit  ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ); // hep-ph/0410050v2
  if ( fSin_Phi_TargPol_LeptonPlane_RF >= 0 )
    fPhi_TargPol_LeptonPlane_RF = fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) );
  if ( fSin_Phi_TargPol_LeptonPlane_RF < 0 )
    fPhi_TargPol_LeptonPlane_RF = 360.0 - std::abs( fRAD2DEG * acos( ( v3QUnitxL.Dot( v3QUnitxS ) ) / ( v3QUnitxL.Mag() * v3QUnitxS.Mag() ) ) );

  return fPhi_TargPol_LeptonPlane_RF;

}

Double_t DEMP_Reaction::Get_Total_Cross_Section(){

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
// 06/09/23 SJDK - Cuts are now ordered as they are applied in the generator
void DEMP_Reaction::Detail_Output(){ // Love Preet changed the order of conditions and added new statements
  
  DEMPDetails << left << setw(70) << "Seed used for the Random Number Generator" << right << setw(20) << fSeed << endl;
  DEMPDetails << endl;
  DEMPDetails << left << setw(70) << "Total events tried" << right << setw(20) << fNGenerated << endl;
  if(UseSolve == true){
    DEMPDetails << left << setw(70) << "Total events cut" << right << setw(20) << (fNaN + fConserve + w_neg_ev + qsq_ev + w_ev  + t_ev + fNSigmaNeg + fNWeightNeg + fSolveEvents_0Sol) << right << setw(20) << ((double) (fNaN + fConserve + w_neg_ev + qsq_ev + w_ev  + t_ev + fNSigmaNeg + fNWeightNeg + fSolveEvents_0Sol)/(double)fNGenerated)*100 << " %" << endl;
    DEMPDetails << left << setw(70) << "Total events recorded" << right << setw(20) << fNRecorded << right << setw(20) << ((double)fNRecorded/(double)fNGenerated)*100 << " %" << endl;
    if (fNGenerated != (fNaN + fConserve + w_neg_ev + qsq_ev + w_ev  + t_ev + fNSigmaNeg + fNWeightNeg + fSolveEvents_0Sol + fNRecorded)){
      DEMPDetails << left << setw(70) << "Total events cut + recorded = events tried?" << right << setw(20) << "NO! ERROR!" << endl;
    }
    else{
      DEMPDetails << left << setw(70) << "Total events cut + recorded = events tried?" << right << setw(20) << "Yes! :)" << endl;
    }
  }
  else{
    DEMPDetails << left << setw(70) << "Total events cut" << right << setw(20) << (fNaN + fConserve + w_neg_ev + qsq_ev + w_ev  + t_ev + fNSigmaNeg + fNWeightNeg) << right << setw(20) << ((double) (fNaN + fConserve + w_neg_ev + qsq_ev + w_ev  + t_ev + fNSigmaNeg + fNWeightNeg)/(double)fNGenerated)*100 << " %" << endl;
    DEMPDetails << left << setw(70) << "Total events recorded" << right << setw(20) << fNRecorded << right << setw(20) << ((double)fNRecorded/(double)fNGenerated)*100 << " %" << endl;
    if (fNGenerated != (fNaN + fConserve + w_neg_ev + qsq_ev + w_ev  + t_ev + fNSigmaNeg + fNWeightNeg + fNRecorded)){
      DEMPDetails << left << setw(70) << "Total events cut + recorded = events tried?" << right << setw(20) << "NO! ERROR!" << endl;
    }
    else{
      DEMPDetails << left << setw(70) << "Total events cut + recorded = events tried?" << right << setw(20) << "Yes! :)" << endl;
    }
  }
  
  DEMPDetails << left << setw(70) << endl << "Cut details -" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to ejectile (X) energy NaN" << right << setw(20) << fNaN << right << setw(20) << ((double)fNaN/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to conservation law check failure" << right << setw(20) << fConserve << right << setw(20) << ((double)fConserve/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to negative Wsq value " << right << setw(20) << w_neg_ev << right << setw(20) << ((double)w_neg_ev/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << Form("Events cut due to qsq < %.1lf or qsq > %.1lf", fQsq_Min, fQsq_Max) << right << setw(20) << qsq_ev << right << setw(20) << ((double)qsq_ev/(double)fNGenerated)*100 << " %" << endl; 
  DEMPDetails << left << setw(70) << Form("Events cut due to W < %.1lf or W > %.1lf", fW_Min, fW_Max) << right << setw(20) << w_ev << right << setw(20) << ((double)w_ev/(double)fNGenerated)*100 << " %" << endl;
  if(UseSolve == true){
    DEMPDetails << left << setw(70) << "Events cut due to solve function finding 0 solutions" << right << setw(20) << fSolveEvents_0Sol << right << setw(20) << ((double)fSolveEvents_0Sol/(double)fNGenerated)*100 << " %" << endl;
  }

  DEMPDetails << left << setw(70) << Form("Events cut due to -t > %.1lf GeV", fT_Max) << right << setw(20) << t_ev << right << setw(20) << ((double)t_ev/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to -ve cross section value" << right << setw(20) << fNSigmaNeg << right << setw(20) << ((double)fNSigmaNeg/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to -ve weight value" << right << setw(20) << fNWeightNeg << right << setw(20) << ((double)fNWeightNeg/(double)fNGenerated)*100 << " %" << endl;
  
  DEMPDetails << left << setw(70) << endl << "Conservation law checks details -" << endl;
  DEMPDetails << left << setw(70) << Form("Total events PASSING conservation law check with tolerance %.2e", fDiff) << right << setw(20) << conserve << endl;
  DEMPDetails << left << setw(70) << "Events cut due to energy conservation check ONLY" << right << setw(20) << ene << right << setw(20) << ((double)ene/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to momentum conservation check ONLY" << right << setw(20) << mom << right << setw(20) << ((double)mom/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to energy AND momentum conservation checks" << right << setw(20) << ene_mom << right << setw(20) << ((double)ene_mom/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to px conservation law check" << right << setw(20) << mom_px << right << setw(20) << ((double)mom_px/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to py conservation law check" << right << setw(20) << mom_py << right << setw(20) << ((double)mom_py/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to pz conservation law check" << right << setw(20) << mom_pz << right << setw(20) << ((double)mom_pz/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to px and py conservation law checks" << right << setw(20) << mom_pxpy << right << setw(20) << ((double)mom_pxpy/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to px and pz conservation law checks" << right << setw(20) << mom_pxpz << right << setw(20) << ((double)mom_pxpz/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to py and pz conservation law checks" << right << setw(20) << mom_pypz << right << setw(20) << ((double)mom_pypz/(double)fNGenerated)*100 << " %" << endl;
  DEMPDetails << left << setw(70) << "Events cut due to px, py and pz conservation law checks" << right << setw(20) << mom_pxpypz << right << setw(20) << ((double)mom_pxpypz/(double)fNGenerated)*100 << " %" << endl;
  
  DEMPDetails << left << setw(70) << endl << "Weight correction factors -" << endl;
  DEMPDetails << left << setw(70) << "Ratio of phase space factors (fPSF_org/fPSF)" << right << setw(20) << ((double)fPSF_org/(double)fPSF) << endl;
  DEMPDetails << left << setw(70) << "Ratio of tried to physical events" << right << setw(20) << ((double)fNGenerated / (double) (fNGenerated - fNaN - fConserve - w_neg_ev))  << endl;
  DEMPDetails << left << setw(70) <<" 1. If the first ratio is not ~1.0 (+/- 0.05), check the number of recorded events. If the number of recorded events is < 50,000, increase the number of attempted events to generate more data. If the the ratio does not converge to ~1.0, multiply all individual event weights by this factor during analysis. Alternatively, this number can be treated as a tolerance on weights."<<endl;
  DEMPDetails << left << setw(70) <<" 2. If the second ratio is large (> XX), then again, multiply all invidiual event weights by this factor during analysis. Again, this could also be considered as a tolerance on weights."<<endl;
   
  DEMPDetails << left << setw(70) << endl << "Energies, angles, and phase space factors -" << endl;
  
  DEMPDetails << right << setw(90) << "User-defined" << right << setw(20) << "Calculated" << right <<  setw(20) << "Actual" << endl;
  
  DEMPDetails << left << setw(70) << "Scattered electron energies (Low)" << right << setw(20) << fScatElec_E_Lo * fElectron_Energy_Col *fm << right << setw(20) << psf_ScatElec_E_min << right << setw(20) << fScatElec_Energy_Col_min * fm << endl;
  
  DEMPDetails << left << setw(70) << "Scattered electron energies (High)" << right << setw(20) << fScatElec_E_Hi * fElectron_Energy_Col *fm << right << setw(20) << psf_ScatElec_E_max << right << setw(20) << fScatElec_Energy_Col_max * fm << endl;
  
  DEMPDetails << left << setw(70) << "Scattered electron angles (Low)" << right << setw(20) << (fScatElec_Theta_I * TMath::RadToDeg()) << right << setw(20) << psf_ScatElec_Theta_min * TMath::RadToDeg() << right << setw(20) << fScatElec_Theta_Col_min * TMath::RadToDeg() << endl; 
  
  DEMPDetails << left << setw(70) << "Scattered electron angles (High)" << right << setw(20) << (fScatElec_Theta_F * TMath::RadToDeg()) << right << setw(20) << psf_ScatElec_Theta_max * TMath::RadToDeg() << right << setw(20) << fScatElec_Theta_Col_max * TMath::RadToDeg() << endl; 
  
  DEMPDetails << left << setw(70) << "Scattered ejectile angles (Low)" << right << setw(20) << (f_Ejectile_Theta_I * TMath::RadToDeg()) << right << setw(20) << psf_Ejectile_Theta_min * TMath::RadToDeg() << right << setw(20) << f_Ejectile_Theta_Col_min * TMath::RadToDeg() << endl;
  
  DEMPDetails << left << setw(70) << "Scattered ejectile angles (High)" << right << setw(20) << (f_Ejectile_Theta_F * TMath::RadToDeg()) << right << setw(20) << psf_Ejectile_Theta_max * TMath::RadToDeg() << right << setw(20) << f_Ejectile_Theta_Col_max * TMath::RadToDeg() << endl;
  
  DEMPDetails << left << setw(70) << "Phase space factors (fPSF fSF_org)" << right << setw(20) << fPSF << right << setw(20) << fPSF_org << endl;
  
  if(UseSolve == true){
    DEMPDetails << left << setw(70) << endl << "Solve function, addtional info -" << endl;
    DEMPDetails << left << setw(70) << "Number of events with 0 Solution" << right << setw(20) << fSolveEvents_0Sol << endl;
    DEMPDetails << left << setw(70) << "Number of events with 1 Solution" << right << setw(20) << fSolveEvents_1Sol << endl;
    DEMPDetails << left << setw(70) << "Number of events with 2 Solution" << right << setw(20) << fSolveEvents_2Sol << endl;
  }
}

////*--------------------------------------------------
/// Functions for different output formats follow

void DEMP_Reaction::Lund_Output(){

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

void DEMP_Reaction::DEMPReact_Pythia6_Out_Init(){

  print_itt = 0;

  DEMPOut << "DEMP Event FILE" << endl;
  DEMPOut << "============================================" << endl;
  DEMPOut << "I, ievent, nParticles, Weight" << endl;
  DEMPOut << "============================================" << endl;
  DEMPOut << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)" << endl;
  DEMPOut << "============================================" << endl;

}

void DEMP_Reaction::DEMPReact_Pythia6_Output(){

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

void DEMP_Reaction::DEMPReact_HEPMC3_Out_Init(){
 
  print_itt = 0;
  DEMPOut << "HepMC::Version 3.02.02" << endl;
  DEMPOut << "HepMC::Asciiv3-START_EVENT_LISTING" << endl;

}

/*--------------------------------------------------*/

void DEMP_Reaction::DEMPReact_HEPMC3_Output(){
  
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
bool DEMP_Reaction::SolnCheck(){

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
double DEMP_Reaction::W_in(){
  return 0;
}

/*--------------------------------------------------*/ 
double DEMP_Reaction::W_out(){
  return 0;
}

/*--------------------------------------------------*/ 

int DEMP_Reaction::Solve(){

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


int DEMP_Reaction::Solve(double theta, double phi){

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

